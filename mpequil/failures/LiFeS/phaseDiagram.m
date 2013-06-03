eps = 1e-10;
A = load('phaseFile.dat');
B = A(:,3:end);
%C = A(:,6:7);

B(find(B<eps)) = NaN;

T = A(:,1)-273.15;
P = A(:,2);
S_s = B(:,1);
S_l = B(:,2);
gas = B(:,3);
FeS2 = B(:,4);


Fe1_xS = B(:,5);
Li3Fe2S4 = B(:,6);
Li2_xFe1_xS2 = B(:,7);
%LiCl_l = C(:,1);
%KCl_l = C(:,2);

figure(1)
plot(1+S_s*0,T,'k*')
hold on
plot(1+S_l*0,T,'b*')
plot(1+gas*0,T,'c*')
plot(2/3+FeS2*0,T,'r*')

%plot(1/2+FeS*0,T,'g*')
if exist('Fe1_xS')
    AAA = load('Fe1_xSPhaseFile.dat');
    BBB = AAA(:,3:end);
    fes = BBB(:,1);
    s = BBB(:,2);
    AAA(find(AAA<eps)) = NaN;
    
    colors = hsv(size(AAA,2)-2);
    figure(3)
    for ( i = 3:size(AAA,2) )
        plot(T,AAA(:,i),'Color',colors(i-2,:),'Marker','*')
        hold on
    end
    hold off 
    figure (1)
    plot((fes+s)./(2*fes+s)+Fe1_xS*0,T,'g*')
else
    plot(1/2+FeS*0,T,'g*')
    plot(8/15+Fe7S8*0,T,'y*')
end
hold off




AA = load('gasPhaseFile.dat');
BB = AA(:,3:end);

AA(find(AA<eps)) = NaN;

colors = hsv(size(AA,2)-2);
figure(2)
for ( i = 3:size(AA,2) )
    plot(T,AA(:,i),'Color',colors(i-2,:),'Marker','*')
    hold on
end
axis([710-273.15 max(T) 0 1]);
legend('S_8','S_7','S_6','S_5','S_4','S_3','S_2','S','location','NorthWest');
hold off


