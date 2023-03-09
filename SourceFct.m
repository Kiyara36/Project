function E = SourceFct(t, InputParas)

if isfield(InputParas,'rep')
    n = floor(t/InputParas.rep); % integer multiples of time that has passed  
    t = t-n*InputParas.rep;
end

if ~ isstruct(InputParas)
    E = InputParas;
else
    E = InputParas.E0*exp(-(t-InputParas.t0)^2/InputParas.wg^2)*exp(1i*(InputParas.we*t + InputParas.phi));
end
end

%t-InputParas.t0 - the difference from t0 from the gaussian wave eqn
%(e^(-x^2/w(width))
%1i (imaginary component for phase shift)
%InputParas.we*t - time which causes modulation
%InputParas.phi - phase shift
%So the equation is giving us a gaussian waveform with a modulation that
%gives us that turning motion. 
% 'rep' - repition time