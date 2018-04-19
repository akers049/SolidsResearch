clear 
clc
%% Transormation matricies for cubic to orthorombic transformation in AuCd

% Parameters for AuCd sources from Bhattacharya Table 4.2 Pg 51
 alpha = 1.0138;
 betta = 0.9491;
 gamma = 1.0350;
 
% Transormation matricies for variants the six varients, sourced from
% Bhattacharya Table 4.2 Pg 51
 
U1 =  [(alpha+gamma)/2,     0,    (alpha-gamma)/2;
               0,          betta,         0      ;
       (alpha-gamma)/2,     0,    (alpha+gamma)/2];

    
U2 =  [(alpha+gamma)/2,     0,    (gamma-alpha)/2;
             0,           betta,         0       ;
      (gamma-alpha)/2       0,    (alpha+gamma)/2];
  
  
U3 =  [(alpha+gamma)/2, (alpha-gamma)/2,    0  ; 
       (alpha-gamma)/2, (alpha+gamma)/2,    0  ;
              0,               0,         betta];
          
U4 =  [(alpha+gamma)/2, (gamma-alpha)/2,    0  ; 
       (gamma-alpha)/2, (alpha+gamma)/2,    0  ;
              0,               0,         betta];
          
U5 =  [betta,         0,               0        ;
         0,    (alpha+gamma)/2,  (alpha-gamma)/2;
         0,    (alpha-gamma)/2,  (alpha+gamma)/2];
    
U6 =  [betta,         0,               0        ;
         0,    (alpha+gamma)/2,  (gamma-alpha)/2;
         0,    (gamma-alpha)/2,  (alpha+gamma)/2];

     
% Put these in a cell array Us containing the transformation matricies, 
% where Us{i} is the transformation matrix for variant i. 

Us = cell(6, 1);
Us = {U1; U2; U3; U4; U5; U6}

save('TransformationMatsForAuCd', 'Us')