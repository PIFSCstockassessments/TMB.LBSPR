#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  // ========================Inputs==================================
  

  // Data in likelihood
  DATA_VECTOR(length_obs);
  DATA_VECTOR(count_obs);
  DATA_VECTOR(GTG_vec);
  DATA_VECTOR(R0_vec);

 
  // Known values
  DATA_SCALAR(M);
  DATA_SCALAR(K);
  DATA_SCALAR(Mat50);
  DATA_SCALAR(Mat95);
  DATA_SCALAR(beta);

  int n_obs_sizeclass = length_obs.size();
  int n_GTG = GTG_vec.size();

  DATA_INTEGER(n_l);

  // =======================Parameters===============================
  //PARAMETER(logF);
  //PARAMETER(logLS50);
  //PARAMETER(logLS95);

  
  PARAMETER(Fmort);
  PARAMETER(LS50);
  PARAMETER(LS95);

  // ====================Transform parameters========================
  //Type F = exp(logF);
  //Type LS50 = exp(logLS50);
  //Type LS95 = exp(logLS95);

    
  // ====================Dynamics=====================================
  Type n_obs=0;
  for(int j=0;j<n_obs_sizeclass;j++){
    n_obs=n_obs+count_obs(j);
  }

  vector<Type> prop_obs(n_obs_sizeclass);
  for(int j=0;j<n_obs_sizeclass;j++){
    prop_obs(j)=count_obs(j)/n_obs;
  }

  vector<Type> L_vec(n_l);
  vector<Type> Sel_vec(n_l);
  vector<Type> M_vec(n_l);
  vector<Type> F_vec(n_l);
  vector<Type> Z_vec(n_l);
  vector<Type> NL_vec(n_l);
  vector<Type> NL_prist_vec(n_l);
  vector<Type> Fec(n_l);
  vector<Type> SSB_vec(n_GTG);
  vector<Type> SSB0_vec(n_GTG);
  vector<Type> DL_vec(n_l);
  matrix<Type> DL_mat(n_l,n_GTG);

  vector<Type> DL_final(n_l);
  vector<Type> DLS_final(n_l);
  vector<Type> prop_expect(n_obs_sizeclass);

  // Populate the length vector
  L_vec(0)=1;
  for(int i=1;i<n_l;i++){ 
    L_vec(i)=L_vec(i-1)+1;  
  }

  // Calculate selectivity at length vector
  for(int i=0;i<n_l;i++){
    Sel_vec(i)=1/(1+exp(-log(19)*(L_vec(i)-LS50)/(LS95-LS50)));
  }

   // Populate M vector
  for(int i=0;i<n_l;i++){
    M_vec(i)=M;

  }  

  // Populate F vector
  for(int i=0;i<n_l;i++){
    F_vec(i)=Fmort*Sel_vec(i);

  }  


  // Populate Z vector
  for(int i=0;i<n_l;i++){
    Z_vec(i)=M_vec(i)+F_vec(i);
  }    


// Start of GTG loop 
for(int GTG=0;GTG<n_GTG;GTG++){  
  
  Type aLinf=GTG_vec(GTG);
  
// Populate Number at length vector
  NL_vec(0)=R0_vec(GTG)*pow((aLinf-L_vec(0)-1)/(aLinf-L_vec(0)),(Z_vec(0)/K));
  for(int i=1;i<n_l;i++){
    NL_vec(i)=NL_vec(i-1)*pow((aLinf-L_vec(i)-1)/(aLinf-L_vec(i)),(Z_vec(i)/K));
    if(L_vec(i)>=(aLinf-1)){NL_vec(i)=0;}
 } 


  // Populate pristine Number at length vector
  NL_prist_vec(0)=R0_vec(GTG)*pow((aLinf-L_vec(0)-1)/(aLinf-L_vec(0)),(M_vec(0)/K));
  for(int i=1;i<n_l;i++){
    NL_prist_vec(i)=NL_prist_vec(i-1)*pow((aLinf-L_vec(i)-1)/(aLinf-L_vec(i)),(M_vec(i)/K));
    if(L_vec(i)>=(aLinf-1)){NL_prist_vec(i)=0;}
  }

  // Calculate fecundity at length vector
  for(int i=0;i<n_l;i++){
    Fec(i)=1/(1+exp(-log(19)*(L_vec(i)-Mat50)/(Mat95-Mat50)))*pow(L_vec(i),beta);
  }

  // Calculate spawning per recruit
  Type aSSB=0;
  for(int i=0;i<n_l-1;i++){
    aSSB=aSSB+1/Z_vec(i)*(NL_vec(i)-NL_vec(i+1))*Fec(i);
  }
  SSB_vec(GTG)=aSSB;

  // Calculate pristine spawning per recruit
  Type aSSB0=0;
  for(int i=0;i<n_l-1;i++){
    aSSB0=aSSB0+1/M_vec(i)*(NL_prist_vec(i)-NL_prist_vec(i+1))*Fec(i);
  }
  SSB0_vec(GTG)=aSSB0;

  // Populate Density at length vector
  for(int i=0;i<n_l-1;i++){
    DL_vec(i)=1/Z_vec(i)*(NL_vec(i)-NL_vec(i+1));
  }

  // Standardize DL_vec column  - THIS WAS A MISTAKE, shouldn't standardize to 1 until after summing DL accross GTGs.
  //Type ColSum=0;
  //for(int i=0;i<n_l;i++){
  //    ColSum=ColSum+DL_vec(i);
  // } 
  
  for(int i=0;i<n_l;i++){
   DL_mat(i,GTG)=DL_vec(i);
   }

} // End of GTG loop

// Sum density by length for all GTGs
for(int i=0;i<n_l;i++){
  DL_final(i)=0;
  for(int GTG=0;GTG<n_GTG;GTG++){
    DL_final(i) = DL_final(i)+DL_mat(i,GTG);
  }
    DL_final(i) = DL_final(i)*Sel_vec(i);
}

  // Standardize DL_final column
  Type ColSum=0;
  for(int i=0;i<n_l;i++){
      ColSum=ColSum+DL_final(i);
  } 
  

  for(int i=0;i<n_l;i++){
   DLS_final(i)=DL_final(i)/ColSum;
   }

// Calculate SPR
Type SSB=0;
Type SSB0=0;   
for(int GTG=0;GTG<n_GTG;GTG++){
   SSB=SSB+SSB_vec(GTG);
   SSB0=SSB0+SSB0_vec(GTG);
}
Type SPR = SSB/SSB0;

 // Re-format to obtain expected number at length
 Type currMinL=0;
 Type currMaxL=length_obs(0);
 int count=0;
 do{
   prop_expect(count)=0;
   for(int i=0;i<n_l;i++){
    if((L_vec(i)>currMinL)&&(L_vec(i)<=currMaxL)){
     prop_expect(count) = prop_expect(count) + DLS_final(i);
        }
   }
   count++;
   currMinL = length_obs(count-1);
   currMaxL = length_obs(count);
   }while(count<n_obs_sizeclass-1);


  Type nll=0;
  Type one_length =0;
  for(int j=0;j<n_obs_sizeclass;j++){
      if(prop_obs(j)>0 && prop_expect(j)>0) one_length=count_obs(j)*log(prop_expect(j)/prop_obs(j));
      //if(prop_obs(j)>0) one_length=count_obs(j)*log(prop_expect(j)/prop_obs(j));
      nll = nll-one_length;
     //Rcout << "Length obs: " << length_obs(j) << " Prop_obs: "<< prop_obs(j) << " Prop_exp: "<< prop_expect(j) << " count_obs: "<< count_obs(j) << std::endl;
  }

  //Rcout << "NLL is " << nll << std::endl;

  REPORT(Fmort);        // plot
  REPORT(LS50);     // plot
  REPORT(LS95);
  REPORT(SPR);
  REPORT(prop_obs);
  REPORT(prop_expect);

  return(nll);
}


