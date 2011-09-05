function result = vadperformance(actual,predicted)

% convert to logical data
actual    = actual==1;
predicted = predicted==1;

% confusion matrix
c = confusionmat(actual,predicted,'order',[0 1]);

% c: 
%    0   1<-predicted
%
%    ?   ?   0 |-> actual
%    ?   ?   1 |

result.speech.recall    = c(2,2)/sum(c(2,:));
result.speech.precision = c(2,2)/sum(c(:,2));

result.nonspeech.recall    = c(1,1)/sum(c(1,:));
result.nonspeech.precision = c(1,1)/sum(c(:,1));

end