function [BiNext]=Bidiag_Francis_Step( Bi )

% Write a single step of a "bidiagonal Francis Step" that introduces the
% bulge and chases it out, working with the bidiagonal matrix in Bi.
% I suggest you do this in stages.
% 
% 1. write it so that you compute the Givens' rotations and you apply
% these to the entire rows and columns to which they are to be applied.  In
% other words, literally do what is in 11.2.4, ignoring where the zeroes
% are when you apply the Givens' rotation.
%
% 2. think carefully about where the nonzeroes are (or are introduced)
% and only update those elements in the rows and/or columns.
%
% Finally, and this is optional, write the routine so you don't corrupt any
% entries below the diagonal and above the first superdiagonal.  Don't
% waste too much time on this.
% Below, we just execute a bunch of these steps, to see if the elements on
% the superdiagonal converge the way they should.

B=Bi;
T=B'*B;
m=size(B,1);
curm=m;
% Set Q

%Q=eye(m,m);
%U=eye(m,m);

T=B'*B;


while curm>1
diff = norm( T(curm,1:curm-1), 1 );
% Compute the 1 norm of the diagonal elements of Ak
diag_1_norm = norm( diag(T), 1 );   
% This stopping criteria needs to be refined.  However, it will be more
% obvious how to do this once we get to the final algorithm
    if diff < 1e-14 * diag_1_norm
% deflate
        curm = curm-1;
        continue
    else
        break
    end
end

if curm==1 % or <3, it should include m==2
   disp("finished, check T")
   BiNext=B;

elseif curm==2
    disp('curm=');
    disp(curm);
% 0. get G0. G0 is from T
    %disp("doing G0,Gh0")
    F=Givens_rotation([T(1,1)-T(curm,curm);T(1,2)]); 
    %G0=eye(m,m);
    %G0(1:2,1:2)=F;
    B(1:2,1:2)=B(1:2,1:2)*F; % A bulge is introduced in B(2,1)   
% 1. Should go right to the last
    F=Givens_rotation([B(1,1);B(2,1)]);
        %Gh=eye(m,m);
        %Gh(i+1:i+2,i+1:i+2)=F;
        %if i+3>m % last Gh
    B(1:2,1:2)=F'*B(1:2,1:2);
    BiNext=B;
    
else
    disp('curm=');
    disp(curm);
% 0. get G0. G0 is from T
    %disp("doing G0,Gh0")
    F=Givens_rotation([T(1,1)-T(curm,curm);T(1,2)]); 
    %G0=eye(m,m);
    %G0(1:2,1:2)=F;
    B(1:2,1:2)=B(1:2,1:2)*F; % A bulge is introduced in B(2,1)   
    
    
% 1. get Ghat0
    F=Givens_rotation([B(1,1);B(2,1)]);
    %Ghat0=eye(m,m);
    %Ghat0(1:2,1:2)=F;
    B(1:2,1:3)=F'*B(1:2,1:3); % A bulge in (1,3)
    %disp('Ghat0*B');
    %disp(B);

%%%% Loop    
    for i=1:curm-2 % from 2 - second to last: m-2 pairs    
        %disp('i=');
        %disp(i);
%2. get Gl
        %disp("doing Gl1,Gh1")
        F=Givens_rotation([B(i,i+1);B(i,i+2)]); % 22, 23?
        %Gl=eye(m,m);
        %Gl(i+1:i+2,i+1:i+2)=F;
        B(i:i+2,i+1:i+2)=B(i:i+2,i+1:i+2)*F; % A bulge should appear in B(3,2)
        %disp('B*Gl');
        %disp(B);
%3. Gh
                   
        F=Givens_rotation([B(i+1,i+1);B(i+2,i+1)]);
        %Gh=eye(m,m);
        %Gh(i+1:i+2,i+1:i+2)=F;
        if i+3>m % last Gh
            B(i+1:i+2,i+1:i+2)=F'*B(i+1:i+2,i+1:i+2);
        else
        B(i+1:i+2,i+1:i+3)=F'*B(i+1:i+2,i+1:i+3);
        %disp('Gh*B');
        %disp(B);
        end
    end
    
    BiNext=B;







end