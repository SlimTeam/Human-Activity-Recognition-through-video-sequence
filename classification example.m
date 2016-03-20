%%%%%%%%%%%%%%%%%%%%% READ FROM EXCEL SHEET
clc;
clear all
% [training]=xlsread('dd.xls'); 
training=xlsread('C:\Users\acer-pc\Desktop\major project 2\training vectors3.xlsx'); 
[r, c]=size(training);
[sample]=xlsread('C:\Users\acer-pc\Desktop\major project 2\testsample1.xlsx'); 
X = training; % use all data for fitting
% Y = species; % response data
Y=[ 1;1;1;1;2;2;2;3;3;3;4;4;4;4;5;5;5;6;6;6] ;
% Y=int2str(Y);
% svmstruct = svmtrain(X,Y);
%  newclass= svmclassify(svmstruct,sample);
%  newclass=int2str(newclass);
class = knnclassify(sample,X, Y, 4,'euclidean');
cp = classperf(Y,class);
cp.CorrectRate 
 disp(class)
% class = classify(sample,X,Y,'diaglinear');
% cp = classperf(testlabel,class);
% cp.CorrectRate 
% disp(class)
% [r, c]=size(training);
% traininglabel = ones(r,1); 
% traininglabel(1:(fix(r/2))) =0; 
% % x=training(:,2:3);
% % x=fix(x);
% % y=traininglabel;
% % training=[{x} {y}];
% % save('111.mat','training')
% [sample]=xlsread('dd1.xls'); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GROUPS
% 
% testlabel=ones(5,1);
% testlabel(1:(fix(5/2))) = 0;
% 
% %%%%%%%%%%%%%%%%%%%%%%% SVM classify
%  T= training(:,(1:2));
%  sample=sample(:,(1:2));
% C=[1;1;1;1;1;0;0;0;0;0;2;2;2;2;2];
% u=unique(C);
% N=length(u);
% if(N>2)
%     disp('multi class problem');
%     itr=1;
%     classes=0;
%     while((classes~=1)&&(itr<=length(u)))
% %this while loop is the multiclass SVM Trick
%     c1=(C==u(itr));
%     newClass=c1;
%     svmStruct = svmtrain(T,newClass,'showplot',true);
%     classes = svmclassify(svmStruct,sample,'showplot',true)
%     itr=itr+1
%     end
% %clc;
% itr=itr-1;
% disp(itr)
% end
 
% svmStruct = svmtrain(T,C,'showplot',true);
% classes = svmclassify(svmStruct,sample,'showplot',true)
% 
% %cp= classperf(traininglabel); 

% 
SVMstruct = svmtrain(X,Y,'kernel_function','polynomial','polyorder',3,'showplot',true); 

crclass=zeros(3);
for indx=1:3
    testname= strcat('C:\Users\acer-pc\Desktop\major project 2\sample vectors\New folder (2)\Book',num2str(indx),'.xlsx');
[sample]=xlsread(testname);
if indx<4
    ch=1;
elseif indx<7
    ch=2;
else
    ch=3;
end

Y1= svmclassify(SVMstruct,sample);
% Z=round(Y);
% res=vec2ind(Z);

res=find(Y1==Y);
crclass(ch,res)=crclass(ch,res)+1;
end
ac=trace(crclass)/3;
disp(strcat('ACCURACY OBTAINED    :' ,num2str(ac),'%'));
beep;

% cp=classperf(Y,newclass);
% % cp.CorrectRate
% 
% idx = (testlabel()==1);
% 
% p = length(testlabel(idx));
% n = length(testlabel(~idx));
% N = p+n;
% 
% tp = sum(testlabel(idx)==class(idx));
% tn = sum(testlabel(~idx)==class(~idx));
% fp = n-tn;
% fn = p-tp;
% 
% tp_rate = tp/p;
% tn_rate = tn/n;
% 
% accuracy = (tp+tn)/N;
% sensitivity = tp_rate;
% specificity = tn_rate;
% precision = tp/(tp+fp);
% recall = sensitivity;
% f_measure = 2*((precision*recall)/(precision + recall));
% gmean = sqrt(tp_rate*tn_rate);
% 
% EVAL = [accuracy sensitivity specificity precision recall f_measure gmean]
%  get(cp)
% disp(newclass)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% DISCRIMINANT classify


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%knn classify
X = training; % use all data for fitting
% Y = species; % response data
 %Y=[1;1;1;1;1;0;0;0;0;0;0;0;0;0;0];
%Y=vertcat('bending','bending','bending','bending','bending','one hand wave','one hand wave','one hand wave','one hand wave','one hand wave','one hand wave','two hand wave','two hand wave','two hand wave','two hand wave','two hand wave');
% Y=int2str(Y);
% mdl = ClassificationKNN.fit(X,Y,'NumNeighbors',4);

%   %cp.countmatrix
%   
% [X,Y] = perfcurve(testlabel,class,2)
% disp(class)




