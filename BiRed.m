function [ A_out, t_out, r_out ] = BiRed( A, t, r )

%% note from the textbook 
%% returns the diagonal and first superdiagonal of the bidiagonal matrix in B, 
%% stores the Householder vectors below the subdiagonal and above the first superdiagonal, 
%% and returns the scalars τ and ρ in vectors t and r
  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );

  [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );
                     
  while ( size( ATL, 1 ) < size( A, 1 ) ) %ends when ATL=A

    [ A00,  a01,        A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,        A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );

    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );

    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );
                                
    %------------------------------------------------------------%
    %disp('size of a12t=')
    %disp(size(a12t))
    %disp('size of A22=')
    %disp(size(A22))
% 1. update alpha11, a21, a12t, A22
% 1.1. do Householder to get rho, tau, and u
    if size(a21)>0 % don't do if there is only alpha11 and no a21
    %   [rho, u2, tau]=Housev(alpha11,a21)
        [u, tau]=Housev1(cat(1,alpha11,a21));
        tau1=u(1);
        u(1)=1;
        %H=eye(size(A22)+1)-u*u'/tau;
     % 1.2 get H
    %H=eye(size(a21,1)+1)-cat(1,1,u2)*cat(1,1,u2)'/tau;
% 1.3 update alpha11, a21
        alpha11=tau1;
        a21=u(2:end);
% 1.4 update (a12t A22)
     %n=size(A22);
     %if n>0
        w12t=(a12t+u(2:end)'*A22)/tau;
        a12t=a12t-w12t;
        A22=A22-u(2:end)*w12t;
     %end
% save rho to t1
        %tau1=rho;
        % do with this (Apr 27)
        tau1=tau;
    
    else
        %alpha11=-alpha11;
        tau1=.5;
  
        
    end
% 2. update A22=A22*H
    if size(a12t,2)>1
% 2.1 get new u12 and rho1
% the textbook uses Housev1 for the row H.
        [u12, rho1]=Housev1(a12t');
% calculate tau
        %disp("tau=")
        %disp(tau)
        tau = ( 1 + u12(2:end)' * u12(2:end) ) / 2;
 
% 2.2 update a12t by storing u12t to a12t
        a12t=u12';
% 2.3 update A22=A22*H(u12, rho1)
% set the first entrey of u12 to one
        u12(1)=1;
        %disp("u12=")
        %disp(u12)
        H=eye(size(A22))-u12*u12'/tau;
        A22=A22*H;
    %else % when there is one superdiagonal only.
     %   a12t=-a12t
     %   rho1=a12t
    
     
    else
        a12t=-a12t;
        rho1=.5;
        
    end  
    
    
    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );

    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );

    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );
 

A_out = [ ATL, ATR
            ABL, ABR ];
  
                                   
                                   
  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  t_out = [ tT
            tB ];
        
  r_out = [ rT
            rB ];

return