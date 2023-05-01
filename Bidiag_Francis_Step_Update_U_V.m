function [ U, B, V ] = Bidiag_Francis_Step_Update_U_V( U_A,B,V_A );

% B: bidiagonal matrix
% U_A, V_A= SVD of B from A
% [ U, B, V ] = Implicit_bidiag_QR_SVD( U, B, V ) calls inside 
%  [ U(:,icur:curm), B(icur:curm,icur:curm), V( :, icur:curm ) ] = 
% Bidiag_Francis_Step_Update_U_V( U(:,icur:curm), B( icur:curm, icur:curm ), V(:,icur:curm ) );

% Now, modify Implicit_bidiag_QR to create Implicit_bidiag_QR_SVD that also
% updates U_A and V_A by applying the Givens' rotations appropriately,
% which means in the end you get the SVD of A, but with the singular values
% not ordered in the correct order.

T=B'*B;
m=size(B,1);
curm=m;

% Set U,V
%U=eye(m,m);
U=U_A;
%V=eye(m,m);
V=V_A;

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

if curm==1
    disp("finished, check T")
    BiNext=B;
    
elseif curm==2
    disp('curm=');
    disp(curm);
% 0. get G0. G0 is from T
    %disp("doing G0,Gh0")
    G0=Givens_rotation([T(1,1)-T(curm,curm);T(1,2)]); 
    %G0=eye(m,m);
    %G0(1:2,1:2)=F;
    B(1:2,1:2)=B(1:2,1:2)*G0; % A bulge is introduced in B(2,1)   
% 1. Should go right to the last
    Ghat0=Givens_rotation([B(1,1);B(2,1)]);
        %Gh=eye(m,m);
        %Gh(i+1:i+2,i+1:i+2)=F;
        %if i+3>m % last Gh
    B(1:2,1:2)=Ghat0'*B(1:2,1:2);
    
    I=eye(m,m);
    I(1:2,1:2)=G0;
    J=eye(m,m);
    J(1:2,1:2)=Ghat0;
    U=U*J; % how to make this faster?
    V=V*I;
    
else
% 0. get G0. G0 is from T
    %disp("doing G0,Gh0")
    G0=Givens_rotation([T(1,1)-T(curm,curm);T(1,2)]); 
    %G0=eye(m,m);
    %G0(1:2,1:2)=F;
    B(1:2,1:2)=B(1:2,1:2)*G0; % A bulge is introduced in B(2,1)   

% 1. get Ghat0
    Ghat0=Givens_rotation([B(1,1);B(2,1)]);
    %Ghat0=eye(m,m);
    %Ghat0(1:2,1:2)=F;
    B(1:2,1:3)=Ghat0'*B(1:2,1:3); % A bulge in (1,3)
    
    I=eye(m,m);
    I(1:2,1:2)=G0;
    J=eye(m,m);
    J(1:2,1:2)=Ghat0;
    U=U*J; % how to make this faster?
    V=V*I;
    %U=U*Ghat0;
    %V=V*G0;

  %%%% Loop    
    for i=1:curm-2 % from 2 - second to last: m-2 pairs    
        %disp('i=');
        %disp(i);
%2. get Gl
        %disp("doing Gl1,Gh1")
        Gl=Givens_rotation([B(i,i+1);B(i,i+2)]); % 22, 23?        
        B(i:i+2,i+1:i+2)=B(i:i+2,i+1:i+2)*Gl; % A bulge should appear in B(3,2)
        
%3. Gh
        Gh=Givens_rotation([B(i+1,i+1);B(i+2,i+1)]);
        if i+3>m % last Gh
            B(i+1:i+2,i+1:i+2)=Gh'*B(i+1:i+2,i+1:i+2);
        else
            B(i+1:i+2,i+1:i+3)=Gh'*B(i+1:i+2,i+1:i+3);
        end
        
        I=eye(m,m);
        I(i+1:i+2,i+1:i+2)=Gl;
        J=eye(m,m);
        J(i+1:i+2,i+1:i+2)=Gh;
        U=U*J; % how to make this faster?
        V=V*I;
        %U=U*Gh;
        %V=V*Gl;

    end  
%end 

end
end