function [precision, recall, cmat] = final_accuracy()

label=ones(43,1);
label(14:24)=2;
label(25:32)=2;
load feat.mat
model=fitcknn(f,label);
out=predict(model,f);

ACTUAL=label;
PREDICTED=out;
idx = (ACTUAL()==1);

p = length(ACTUAL(idx));
n = length(ACTUAL(~idx));
N = p+n;

tp = sum(ACTUAL(idx)==PREDICTED(idx));
tn = sum(ACTUAL(~idx)==PREDICTED(~idx));
fp = n-tn;
fn = p-tp;

tp_rate = tp/p;
tn_rate = tn/n;

accuracy = (tp+tn)/N;
sensitivity = tp_rate;
specificity = tn_rate;
precision = tp/(tp+fp);
recall = sensitivity;
f_measure = 2*((precision*recall)/(precision + recall));
gmean = sqrt(tp_rate*tn_rate);

EVAL = [accuracy sensitivity specificity precision recall f_measure gmean];
fprintf('\n accuracy  classifer = %f\n',accuracy);
fprintf('\n sensitivity  classifer = %f\n',sensitivity);
fprintf('\n specificity  classifer = %f\n',specificity);
fprintf('\n precision classifer = %f\n',precision);
fprintf('\n recall  classifer = %f\n',recall);
end
