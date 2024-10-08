function [indout,offset] = hPRACHDetect(carrier,prach,waveform,indin) 
%hPRACHDetect Physical random access channel (PRACH) detection
%   [INDOUT,OFFSET] = hPRACHDetect(CARRIER,PRACH,WAVEFORM,INDIN) performs
%   PRACH detection given carrier-specific configuration CARRIER,
%   PRACH-specific configuration PRACH, received signal potentially
%   containing a PRACH transmission WAVEFORM, and range of preamble indices
%   for which to search specified in INDIN. The detector performs each
%   distinct correlation required to cover all preamble indices specified
%   in INDIN, and searches the output of the correlations for peaks which
%   exceed a detection threshold. The position of the peak in the
%   correlator output is used to determine the preamble index that was
%   detected and its associated timing offset, with the preamble index and
%   timing offset being returned in INDOUT and OFFSET respectively.
%
%   CARRIER is a carrier-specific configuration object, as described in
%   <a href="matlab:help('nrCarrierConfig')" >nrCarrierConfig</a> with these properties:
%
%   SubcarrierSpacing   - Subcarrier spacing in kHz (15, 30, 60, 120, 240)
%   CyclicPrefix        - Cyclic prefix ('normal', 'extended')
%   NSizeGrid           - Number of resource blocks in carrier resource
%                         grid (1...275)
%
%   PRACH is a PRACH-specific configuration object, as described in
%   <a href="matlab:help('nrPRACHConfig')" >nrPRACHConfig</a> with these properties:
%
%   FrequencyRange      - Frequency range (used in combination with
%                         DuplexMode to select a configuration table from
%                         TS 38.211 Table 6.3.3.2-2 to 6.3.3.2-4)
%   DuplexMode          - Duplex mode (used in combination with
%                         FrequencyRange to select a configuration table
%                         from TS 38.211 Table 6.3.3.2-2 to 6.3.3.2-4)
%   ConfigurationIndex  - Configuration index, as defined in TS 38.211
%                         Tables 6.3.3.2-2 to 6.3.3.2-4
%   SubcarrierSpacing   - PRACH subcarrier spacing in kHz
%   NPRACHSlot          - PRACH slot number
%
%   WAVEFORM is an N-by-P matrix containing the received time-domain signal
%   in which to search for PRACH transmissions, where N is the number of
%   time domain samples and P is the number of receive antennas. Such a
%   waveform is generated by the <matlab:help('hNRPRACHWaveformGenerator.m')
%   hNRPRACHWaveformGenerator.m> helper function, with P=1.
%   It is assumed that any PRACH signal in WAVEFORM is synchronized such
%   that the first sample of WAVEFORM corresponds to the start of an uplink
%   subframe, therefore any delay from the start of WAVEFORM to the first
%   sample of the PRACH therein will be interpreted by the detector as a
%   timing offset.
%
%   INDIN is a vector of preamble indices within the cell for which to
%   search. INDIN can be between 1 and 64 in length, containing values
%   between 0 and 63.
%
%   Output INDOUT is the index (included in INDIN) indicating which
%   preamble was found during the search. If no index resulted in a
%   correlation above the detection threshold or the maximum correlation
%   was obtained for an index not included in INDIN, INDOUT is empty.
%
%   Output OFFSET is the timing offset expressed in samples at the input
%   sampling rate. The timing offset estimate has an integer part
%   corresponding to the correlation peak sample position and a fractional
%   part estimating the fractional delay present in the correlation peak.
%   The cyclic shift in the frequency domain present in the PRACH preamble
%   can contribute to this fractional delay. If no index from INDIN results
%   in a correlation above the detection threshold or the maximum
%   correlation was obtained for an index not included in INDIN, OFFSET is
%   empty.
%
%   The detector first calls [~,INFO]=nrPRACH(CARRIER,PRACH) to establish
%   the set of root sequences INFO.RootSequence required to cover all
%   preamble indices in INDIN. A correlation is then performed for each
%   distinct value in INFO.RootSequence, with the inputs to the correlation
%   being the input WAVEFORM and a locally generated PRACH waveform. The
%   correlation is performed in the frequency domain - multiplication of
%   the FFT of the useful part of the locally generated PRACH waveform by a
%   portion of the input WAVEFORM extracted with the same timing as the
%   useful part of the locally generated PRACH waveform, followed by an
%   IFFT to give the correlation. Further fields from INFO are then used to
%   establish the length of the window of the correlator output that
%   corresponds to each preamble index (the zero correlation zone). The
%   preamble index is established by testing the position of the peak in
%   the correlator output to determine if it lies in the window of the
%   correlator output given by the cyclic shift for each preamble index in
%   turn, with the offset within the window being used to compute the
%   timing offset.
%
%   Example:
%   % Detect a PRACH preamble format 0 which has been delayed by 7 samples.
%   % Note that the timing offset includes a fractional part which is an 
%   % estimate of the fractional delay present in the correlation peak. 
%   % This is due to the cyclic shift present in the PRACH preamble. A 
%   % cyclic shift in the frequency domain is a delay in the time domain.
%   
%   carrier = nrCarrierConfig;
%   carrier.NSizeGrid = 6;
%   prach = nrPRACHConfig;
%   prach.ConfigurationIndex = 27;
%   prach.ZeroCorrelationZone = 1;
%   prach.PreambleIndex = 44;
%   waveconfig.NumSubframes = 1;
%   waveconfig.Carriers = carrier;
%   waveconfig.PRACH.Config = prach;
%   tx = hNRPRACHWaveformGenerator(waveconfig);
%   rx = [zeros(7,1); tx]; % delay PRACH 
%   [index,offset] = hPRACHDetect(carrier,prach,rx,(0:63).')
%
%   See also nrPRACHConfig, nrPRACH, nrPRACHIndices, nrPRACHGrid,
%   nrPRACHOFDMInfo, nrPRACHOFDMModulate, hNRPRACHWaveformGenerator.

%   Copyright 2019-2022 The MathWorks, Inc.
    
    % Pre-configure empty outputs
    indout = [];
    offset = [];
    
    % Store most used PRACH properties for ease of access
    LRA = prach.LRA;
    subcarrierSpacing = prach.SubcarrierSpacing*1e3; % Subcarrier spacing in Hz
    prachDuration = prach.PRACHDuration;
    
    % For format C2, the last OFDM symbol in each time occasion is empty,
    % and therefore the length of the preamble is one symbol shorter than 
    % the value of PRACHDuration
    if strcmpi(prach.Format,'C2')
        prachDuration = prachDuration - 1;
    end
    
    % Extract PRACH-related information |prachInfo| and OFDM-related
    % information for PRACH |prachOFDMInfo|
    [~,prachInfo] = arrayfun(@(x)nrPRACH(carrier,setfield(prach,'PreambleIndex',x)),0:63,'UniformOutput',false);
    prachInfo = [prachInfo{:}]; % Convert prachInfo to an array of struct
    ofdmInfo = nrPRACHOFDMInfo(carrier,prach,'Windowing',0);
    
    % Find set of root sequences required for set of input preamble indices
    u = unique([prachInfo(indin+1).RootSequence]);
    
    % Set up of number of receive antennas |NRxAnts|, correlation outputs
    % |c| first preamble index for each distinct root sequence
    % |preambleidx| and detection threshold for each antenna |threshold|.
    NRxAnts = size(waveform,2);
    c = cell(length(u),1);
    preambleidx = zeros(length(u),1);
    threshold = zeros(NRxAnts,1);
    
    % Configure correlator input dimensions
    numOFDMSymbPerSlot = numel(ofdmInfo.SymbolLengths);
    symbLoc = prach.SymbolLocation - numOFDMSymbPerSlot*prach.ActivePRACHSlot;
    start = ofdmInfo.OffsetLength + sum(ofdmInfo.SymbolLengths(1:symbLoc));
    if (numOFDMSymbPerSlot > 0)
        start = start + ofdmInfo.CyclicPrefixLengths(symbLoc+1);
        duration = ofdmInfo.Nfft;
    else
        duration = 0;
    end
    ncorrs = prachDuration;
    
    % Store the configuration parameters needed to generate the reference PRACH waveform
    waveconfig.NumSubframes = prach.SubframesPerPRACHSlot;
    waveconfig.Windowing = 0;
    waveconfig.Carriers = carrier;
    waveconfig.PRACH.Config = prach;
    
    % Ratio of carrier and PRACH subcarrier spacings
    % TS 38.211 Section 5.3.2
    K = carrier.SubcarrierSpacing  / prach.SubcarrierSpacing;
    
    % Perform all necessary correlations i.e. those with distinct root sequences
    for i = 1:length(u)

        % Find first preamble index preambleidx(i) to give root sequence u(i).
        preambleidx(i) = find([prachInfo.RootSequence]==u(i),1)-1;
        waveconfig.PRACH.Config.PreambleIndex = preambleidx(i);
        info = prachInfo(preambleidx(i)+1);
        
        % Generate reference PRACH sequence
        refPRACH = hNRPRACHWaveformGenerator(waveconfig);
        refPRACH = refPRACH(start+(1:duration));
        
        % For each receive antenna
        for p = 1:NRxAnts
            
            % Perform correlation(s) of input on the pth antenna with
            % reference and store result
            cp = zeros(duration,1);
            for k = 1:ncorrs
                rx = waveform(start+((k-1)*duration)+(1:duration),p);
                rxFFT = fft(rx);                        
                cp = cp+abs(ifft(rxFFT.*conj(fft(refPRACH)))).^2;
            end
            cp = cp/sqrt(ncorrs);
            
            % Form threshold for subsequent detection based on estimation
            % of the received signal power on pth antenna within PRACH
            % bandwidth. The ratio K/12 deals with difference between
            % expected level of correlation output for different PRACH and
            % carrier subcarrier spacings. The scaling factor of 166 has
            % been determined empirically to meet the probability of
            % detection in TDL-C.
            if (i==1)
                mask = abs(fft(refPRACH))>0.1;
                threshold(p) = var(ifft(rxFFT.*mask)) * 100 * (K/12);
            end
            
            % Combine correlations from each antenna
            if (p==1)
                c{i} = cp;
            else
                c{i} = c{i}+cp;
            end            
            
        end
        c{i} = c{i}/NRxAnts;
        
        % For restricted sets, determine the size of the cyclic offset due
        % to Doppler shift, in order to deal with the side peaks caused by
        % loss of orthogonality due to high Doppler shift. The side peaks
        % are dealt with by combining each interval of length
        % |cyclicOffset| in the correlation output.
        if ~strcmpi(prach.RestrictedSet,'UnrestrictedSet')
            cyclicOffset=fix((info.CyclicOffset/LRA)*ofdmInfo.SampleRate/subcarrierSpacing);
            x=c{i};
            c{i}=(x((1+cyclicOffset):end)+x(1:(end-cyclicOffset)))/sqrt(2);
        end
        
    end
    
    % Combine threshold across antennas
    threshold = mean(threshold)/sqrt(NRxAnts);
    
    % Determine the length of the zero correlation zone
    zcz = (prachInfo(1).NumCyclicShifts/LRA)*ofdmInfo.SampleRate/subcarrierSpacing;
    
    % The following parameter specifies the fraction of the timing window
    % at the end of the timing window for one preamble that will be
    % considered as belonging to the next preamble and having a timing
    % offset of zero. This effectively excludes timing offsets of above
    % (1.0-deadzone) of the maximum and ensures detection of preambles with
    % low timing offset where noise has caused the peak of the correlation
    % to be slightly into the previous preamble's timing window. The value
    % configured below corresponds to the duration of the main lobe of the
    % autocorrelation of the PRACH. Zero is used for the case in which NCS
    % is 0 as there is only one preamble per correlation.
    if (zcz~=0)
        deadzone = ofdmInfo.SampleRate/(LRA*subcarrierSpacing)/zcz;
    else
        deadzone = 0;
    end
    
    % Detect preambles. The implementation here will detect a single
    % preamble index with the strongest correlation across all correlators,
    % provided a detection threshold is exceeded. 
    linearidx = 0:duration-1;
    bestCorr = 0;
    
    % For each unique root sequence
    for i = 1:length(u)
        
        % If the maximum value of the correlation for this root sequence
        % exceeds the detection threshold and is the highest maximum value
        % across correlations checked so far
        if (max(c{i})>threshold)
            if (max(c{i})>bestCorr)
                
                % Record the peak correlation value and its position
                bestCorr = max(c{i});
                maxpos = mod(linearidx(c{i}==bestCorr)+(deadzone*zcz),length(c{i}))-(deadzone*zcz);
                
                % Establish the preamble index and timing offset from the
                % correlation peak position
                if (info.NumCyclicShifts==0)
                    indout = preambleidx(i);
                    offset = maxpos;
                else
                    % Find the set of cyclic shifts v = 0...maxv for this
                    % root sequence.
                    maxv = find([prachInfo.RootSequence]==u(i),1,'Last')-preambleidx(i)-1;
                    cyclicShift = (mod(LRA-[prachInfo(preambleidx(i)+(1:(maxv+1))).CyclicShift],LRA)/LRA)*ofdmInfo.SampleRate/subcarrierSpacing;
                    
                    % Establish the set of offsets from the peak
                    % correlation position to the set of cyclic shifts for
                    % this root sequence
                    offsetpos = maxpos-cyclicShift;
                    if ~strcmpi(prach.RestrictedSet,'UnrestrictedSet')
                        % For restricted sets, determine the size of the
                        % cyclic offset due to Doppler shift, in order to
                        % deal with the side peaks caused by loss of
                        % orthogonality due to high Doppler shift
                        cyclicOffset=(info.CyclicOffset/LRA)*ofdmInfo.SampleRate/subcarrierSpacing;
                        offsetpos=[offsetpos-cyclicOffset offsetpos offsetpos+cyclicOffset]; %#ok<AGROW>
                    end
                    
                    % Find the value of v for the detected preamble and
                    % compute the final preamble index
                    vdash = find(floor(offsetpos/zcz + deadzone)==0)-1;
                    if (~isempty(vdash))
                        vdash = vdash(1);
                    end
                    v = mod(vdash,maxv+1);
                    indout = preambleidx(i)+v;
                    
                    % Establish the timing offset from the correlation peak position
                    offset = offsetpos(vdash+1);
                    offset = max(offset, 0);
                end
            end
        end
    end
    
    % Remove output values if detected preamble is not part of input
    % preamble index range INDIN
    if(~isempty(indout))
        if(isempty(find(indin==indout,1)))
            indout = [];
            offset = [];
        end
    end
end