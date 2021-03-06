function y = sigmoidCurveModel(x,a,b,c,d)
% FITSIGMOIDCURVE fits a generalized sigmoid function 
% to the input data: a = lower asymptote, b = upper asymptote
% c = usually 1, d = growth rate, k = linked to mid value
y = zeros(size(x));

% for i = 1:length(x)
    y=a + (b-a) ./ (1 + 10.^((c-x)*d));
% end
end