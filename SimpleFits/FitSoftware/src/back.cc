#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include <iostream>


DiTauConstrainedFitter::DiTauConstrainedFitter(LorentzVectorParticle TauA1,TrackParticle MuTrack, TVector3 PVertex, TMatrixTSym<double> VertexCov):
  LagrangeMultipliersFitter(),
  Zwidth_(2.5*2.5),
  Zmass_(91.5)
{
//        TRandom *ran = new TRandom();  //  Generate Random value around z mass
//        double ZMs = ran->Gaus(91.18,2.5);

  //   std::cout<<"RandomlyGeneratedZmass1 _"<<RandomlyGeneratedZmass_<<"  ZMs " <<ZMs <<std::endl;

  // LorentzVectorParticle Mu = ConvertTrackParticleToLorentzVectorParticle(MuTrack);
  LorentzVectorParticle TauMuGuess  = EstimateTauMu( PVertex, VertexCov, TauA1.Vertex(),TauA1.VertexCov(), MuTrack,  TauA1,Zmass_ );


  particles_.push_back(TauA1);
  particles_.push_back(TauMuGuess);

//   std::cout<<"DiTauConstrainedFitter  TauMuGuess"<<   TauA1.LV().Px() <<"  "<<TauA1.LV().Py()<<"  "<<TauA1.LV().Pz()<< std::endl;
   std::cout<<"DiTauConstrainedFitter  TauMuGuess"<<   TauMuGuess.LV().Px() <<"  "<<TauMuGuess.LV().Py()<<"  "<<TauMuGuess.LV().Pz()<< std::endl;

//  std::cout<<"checkMass2  "<<(TauMuGuess.LV()+TauA1.LV()).M()<<std::endl;


  isconfigured=false;

  // setup 6 by 6 matrix 
  //  int size=particles_.size()*3 + 1;   //  3*momenta  + z mass

  int size=particles_.size()*3;   //  3*momenta  + z mass
  TMatrixT<double>    inpar(size,1);
  TMatrixTSym<double> incov(size);
 
  // Get primary vertex information
  if(VertexCov.GetNrows()!=LorentzVectorParticle::NVertex)return;


  //set input paramterts:  TauA1 - TauMu
   inpar(taua1_px,0)=TauA1.LV().Px();
   inpar(taua1_py,0)=TauA1.LV().Py();
   inpar(taua1_pz,0)=TauA1.LV().Pz();

   inpar(taumu_px,0)=TauMuGuess.LV().Px();
   inpar(taumu_py,0)=TauMuGuess.LV().Py();
   inpar(taumu_pz,0)=TauMuGuess.LV().Pz();

   //   inpar(z_m,0)=Zmass_;

   std::cout<<"mass check2   "<< (TauA1.LV()+TauMuGuess.LV()).M() <<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<<   inpar(LorentzVectorParticle::vx,0)<<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<<   inpar(LorentzVectorParticle::vy,0)<<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<<   inpar(LorentzVectorParticle::vz,0)<<std::endl;
   
   std::cout<<"LorentzVectorParticle 2 "<<   inpar(LorentzVectorParticle::px,0)<<std::endl;
   std::cout<<"LorentzVectorParticle 2 "<<   inpar(LorentzVectorParticle::py,0)<<std::endl;
   std::cout<<"LorentzVectorParticle 2 "<<   inpar(LorentzVectorParticle::pz,0)<<std::endl;

   std::cout<<"LorentzVectorParticle 1 "<< TauA1.LV().Px()  <<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<< TauA1.LV().Py()   <<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<< TauA1.LV().Pz()   <<std::endl;
   
   std::cout<<"LorentzVectorParticle 2 "<< TauMuGuess.LV().Px()   <<std::endl;
   std::cout<<"LorentzVectorParticle 2 "<< TauMuGuess.LV().Py()   <<std::endl;
   std::cout<<"LorentzVectorParticle 2 "<< TauMuGuess.LV().Pz()   <<std::endl;


    int TauA1offset=0;
    int TauMuoffset=3;
   
    for(int i=0; i<3;i++){
     for(int j=0; j<3;j++){
       incov(i+TauA1offset,j+TauA1offset)=TauA1.Covariance(i+3,j+3); 
       incov(i+TauMuoffset,j+TauMuoffset)=TauMuGuess.Covariance(i+3,j+3);
     } 
    }
//     incov(0,6) = 0;       incov(6,0) = 0; 
//     incov(1,6) = 0;       incov(6,1) = 0; 
//     incov(2,6) = 0;       incov(6,2) = 0; 
//     incov(3,6) = 0;       incov(6,3) = 0; 
//     incov(4,6) = 0;       incov(6,4) = 0; 
//     incov(5,6) = 0;       incov(6,5) = 0; 

//    incov(6,6) = Zwidth_;
    std::cout<<"inpar "<<std::endl;
   for(int str =0; str < inpar.GetNrows(); str++){
     for(int kol =0; kol < inpar.GetNcols(); kol++){
        std::cout<<"  "<< inpar(str,kol)<<"  ";

     }    
      std::cout<<std::endl;
   }

    std::cout<<"TauMuGuessCovariance "<<std::endl;
    for(int str =0; str < TauMuGuess.VertexCov().GetNrows(); str++){
     for(int kol =0; kol < TauMuGuess.VertexCov().GetNcols(); kol++){
        std::cout<<"  "<< TauMuGuess.VertexCov()(str,kol)<<"  ";

     }    
      std::cout<<std::endl;
   }




 // store expanded par for computation of final par (assumes fit has neglegible impact on a1 correlations with vertex uncertainties)


  //----
//    std::cout<<"DiTauConstraint inpar "<<std::endl;
//   for(int str =0; str < inpar.GetNrows(); str++){
//     for(int kol =0; kol < inpar.GetNcols(); kol++){
//        std::cout<<"  "<< inpar(str,kol)<<"  ";

//     }    
//      std::cout<<std::endl;
//   }

//   //----
//     std::cout<<"DiTauConstraint incov "<<std::endl;
//    for(int str =0; str < incov.GetNrows(); str++){
//      for(int kol =0; kol < incov.GetNcols(); kol++){
//         std::cout<<"  "<< incov(str,kol)<<"  ";

//      }    
//       std::cout<<std::endl;
//    }


   exppar.ResizeTo(size,1);
   exppar=ComputeInitalExpPar(inpar);
   expcov.ResizeTo(size,size);
   expcov=ErrorMatrixPropagator::PropogateError(&DiTauConstrainedFitter::ComputeInitalExpPar,inpar,incov);
   //

//    exppar = inpar;
//    expcov = incov;
   ////----
//     std::cout<<"DiTauConstraint exppar "<<std::endl;
//    for(int str =0; str < exppar.GetNrows(); str++){
//      for(int kol =0; kol < exppar.GetNcols(); kol++){
//         std::cout<<"  "<< exppar(str,kol)<<"  ";

//      }    
//      std::cout<<std::endl;
//    }

// //     //----
//     std::cout<<"DiTauConstraint expcov "<<std::endl;
//    for(int str =0; str < expcov.GetNrows(); str++){
//      for(int kol =0; kol < expcov.GetNcols(); kol++){
//         std::cout<<"  "<< expcov(str,kol)<<"  ";

//      }    
//       std::cout<<std::endl;
//    } 



  // store linearization point
  TMatrixT<double> PAR_0(size,1);
  par_0.ResizeTo(size);
  cov_0.ResizeTo(size,size);
  PAR_0=ComputeExpParToPar(exppar);
  for(int i=0; i<npar;i++)par_0(i)=PAR_0(i,0);
  cov_0=ErrorMatrixPropagator::PropogateError(&DiTauConstrainedFitter::ComputeExpParToPar,exppar,expcov);

  for(int i=0; i<npar;i++){
    for(int j=0;j<npar;j++){cov_0(i,j)=expcov(i,j);}
  }
  // set up inital point for fit (cov handled in Fit() function)
  par.ResizeTo(npar);
  par=par_0;

// std::cout<<"DiTauConstraint  cov_0"<<std::endl;
// for(int str =0; str < cov_0.GetNrows(); str++){
//   for(int kol =0; kol < cov_0.GetNcols(); kol++){
//     std::cout<<"  "<< cov_0(str,kol)<<"  ";
    
//   }    
//   std::cout<<std::endl;
//  } 



// std::cout<<"DiTauConstraint  par_0"<<std::endl;
// for(int str =0; str < par.GetNrows(); str++){
//   for(int kol =0; kol < par.GetNcols(); kol++){
//     std::cout<<"  "<< par(str,kol)<<"  ";
    
//   }    
//   std::cout<<std::endl;
//  } 



//   TLorentzVector a1(par(a1_px),par(a1_py),par(a1_pz),sqrt(par(a1_m)*par(a1_m)+par(a1_px)*par(a1_px)+par(a1_py)*par(a1_py)+par(a1_pz)*par(a1_pz)));
//   double phi(par(tau_phi)),theta(par(tau_theta));

  isconfigured=true; 
//   std::cout<<"deb 2"<<std::endl; 
}


TMatrixT<double> DiTauConstrainedFitter::ComputeInitalExpPar(TMatrixT<double> &inpar){
  // std::cout<<"deb ComputeInitalExpPar "<<std::endl; 
  TMatrixT<double> outpar(6,1);
  int offset=0;//LorentzVectorParticle::NVertex;// for TauA1
  outpar(taua1_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taua1_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taua1_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  offset+=3;//LorentzVectorParticle::NLorentzandVertexPar; // for TauNu

  outpar(taumu_px,0)=inpar(LorentzVectorParticle::vx+offset,0);
  outpar(taumu_py,0)=inpar(LorentzVectorParticle::vy+offset,0);
  outpar(taumu_pz,0)=inpar(LorentzVectorParticle::vz+offset,0);

  // outpar(z_m,0)=inpar(6,0);

  //  std::cout<<"outpar(taua1_px,0)  "<<outpar(taua1_px,0)<<"  "<<inpar(LorentzVectorParticle::vx+offset,0)  <<std::endl; 
  return outpar; 

}


TMatrixT<double> DiTauConstrainedFitter::ComputeExpParToPar(TMatrixT<double> &inpar){
//  std::cout<<"deb 5"<<std::endl; 
  TMatrixT<double> outpar(npar,1);
  for(int i=0;i<npar;i++){outpar(i,0)=inpar(i,0);}
//  std::cout<<"deb 6"<<std::endl; 
  return outpar;
}


TMatrixT<double> DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);

   outpar(3,0)=inpar(3,0);
   outpar(4,0)=inpar(4,0);
   outpar(5,0)=inpar(5,0);
   outpar(6,0)=1.777;

//    std::cout<<"--ComputeTauMuLorentzVectorPar  "<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(0,0)"<< outpar(3,0)<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(1,0)"<< outpar(4,0)<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(2,0)"<< outpar(5,0)<<std::endl; 
//    std::cout<<"--ComputeTauMuLorentzVectorPar  outpar(3,0)"<< outpar(6,0)<<std::endl; 


  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar(TMatrixT<double> &inpar){
  //  start with index 3 to fill only momenta part
  TMatrixT<double> outpar(7,1);
  outpar(3,0)=inpar(0,0);
  outpar(4,0)=inpar(1,0);
  outpar(5,0)=inpar(2,0);
  outpar(6,0) =1.777;
//    std::cout<<"ComputeTauA1LorentzVectorPar  "<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(0,0)"<< outpar(3,0)<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(1,0)"<< outpar(4,0)<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(2,0)"<< outpar(5,0)<<std::endl; 
//    std::cout<<"ComputeTauA1LorentzVectorPar  outpar(3,0)"<< outpar(6,0)<<std::endl; 
  return outpar;
}

TMatrixT<double> DiTauConstrainedFitter::ComputeMotherLorentzVectorPar(TMatrixT<double> &inpar){
//   std::cout<<"deb 11"<<std::endl;     
  TMatrixT<double> outpar(7,1);
  TMatrixT<double> Taunupar=ComputeTauMuLorentzVectorPar(inpar);
  TMatrixT<double> Taua1par=ComputeTauA1LorentzVectorPar(inpar);

  outpar(3,0)=Taunupar(3,0)+Taua1par(3,0);
  outpar(4,0)=Taunupar(4,0)+Taua1par(4,0);
  outpar(5,0)=Taunupar(5,0)+Taua1par(5,0);

  double Etaumu2=pow(Taunupar(3,0),2.0)+pow(Taunupar(4,0),2.0)+pow(Taunupar(5,0),2.0)+pow(Taunupar(6,0),2.0);
  double Etaua12=pow(Taua1par(3,0),2.0)+pow(Taua1par(4,0),2.0)+pow(Taua1par(5,0),2.0)+pow(Taua1par(6,0),2.0);
  double P2=pow(outpar(3,0),2.0)+pow(outpar(4,0),2.0)+pow(outpar(5,0),2.0);
  outpar(6,0)=sqrt(fabs(pow(sqrt(Etaumu2)+sqrt(Etaua12),2.0)-P2));
//   std::cout<<"outpar(0,0) "<<outpar(3,0) <<std::endl;     
//   std::cout<<"outpar(1,0) "<<outpar(4,0) <<std::endl;     
//   std::cout<<"outpar(2,0) "<<outpar(5,0) <<std::endl;     
//   std::cout<<"outpar(3,0) "<<outpar(6,0) <<std::endl;     
  return outpar;
}

void DiTauConstrainedFitter::UpdateExpandedPar(){
  // assumes changes to a1 correlation to vertex is small
  if(par.GetNrows()==npar && cov.GetNrows() && exppar.GetNrows()==npar && expcov.GetNrows()) return;
  for(int i=0; i<npar;i++){
    exppar(i,0)=par(i);
    for(int j=0; j<npar;j++){expcov(i,j)=cov(i,j);}
  }
}



std::vector<LorentzVectorParticle> DiTauConstrainedFitter::GetReFitDaughters(){
  std::vector<LorentzVectorParticle> refitParticles;
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> a1=ComputeTauA1LorentzVectorPar(exppar);

//    std::cout<<"a1 cross ComputeTauA1LorentzVectorPar  "<<std::endl; 
//    std::cout<<"a1 cross ComputeTauA1LorentzVectorPar  outpar(0,0)"<< a1(3,0)<<std::endl; 
//    std::cout<<"a1 cross ComputeTauA1LorentzVectorPar  outpar(1,0)"<< a1(4,0)<<std::endl; 
//    std::cout<<"a1 cross ComputeTauA1LorentzVectorPar  outpar(2,0)"<< a1(5,0)<<std::endl; 
//    std::cout<<"a1 cross ComputeTauA1LorentzVectorPar  outpar(3,0)"<< a1(6,0)<<std::endl; 

  TMatrixTSym<double> a1cov=ErrorMatrixPropagator::PropogateError(&DiTauConstrainedFitter::ComputeTauA1LorentzVectorPar,exppar,expcov);
  refitParticles.push_back(LorentzVectorParticle(a1,a1cov,fabs(PDGInfo::tau_plus)*c,c,b));

//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  "<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(0,0)"<< refitParticles.at(0).LV().Px()<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(1,0)"<< refitParticles.at(0).LV().Py()<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(2,0)"<< refitParticles.at(0).LV().Pz()<<std::endl; 
//     std::cout<<"a1 cross2 ComputeTauA1LorentzVectorPar  outpar(3,0)"<< refitParticles.at(0).LV().M()<<std::endl; 

  TMatrixT<double> nu=ComputeTauMuLorentzVectorPar(exppar);
  nu(0,0)= particles_.at(1).Parameter(LorentzVectorParticle::vx);
  nu(1,0)= particles_.at(1).Parameter(LorentzVectorParticle::vy);
  nu(2,0)= particles_.at(1).Parameter(LorentzVectorParticle::vz);
  TMatrixTSym<double> nucov=ErrorMatrixPropagator::PropogateError(&DiTauConstrainedFitter::ComputeTauMuLorentzVectorPar,exppar,expcov);
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      nucov(i,j)=particles_.at(1).VertexCov()(i,j);
      
    }
  }
  
  refitParticles.push_back(LorentzVectorParticle(nu,nucov,PDGInfo::tau_minus,0.0,b));
//    std::cout<<"a1 cross3 ComputeTauA1LorentzVectorPar  "<<std::endl; 
//    std::cout<<"a1 cross3 ComputeTauA1LorentzVectorPar outpar(0,0)"<< refitParticles.at(1).LV().Px()<<std::endl; 
//    std::cout<<"a1 cross3 ComputeTauA1LorentzVectorPar  outpar(1,0)"<< refitParticles.at(1).LV().Py()<<std::endl; 
//    std::cout<<"a1 cross3 ComputeTauA1LorentzVectorPar  outpar(2,0)"<< refitParticles.at(1).LV().Pz()<<std::endl; 
//    std::cout<<"a1 cross3 ComputeTauA1LorentzVectorPar  outpar(3,0)"<< refitParticles.at(1).LV().M()<<std::endl; 



  return refitParticles; 
}

LorentzVectorParticle DiTauConstrainedFitter::GetMother(){
  UpdateExpandedPar();
  double c(0),b(0);
  for(unsigned int i=0;i<particles_.size();i++){c+=particles_.at(i).Charge();b=particles_.at(i).BField();}
  TMatrixT<double> m=ComputeMotherLorentzVectorPar(exppar);
  TMatrixTSym<double> mcov=ErrorMatrixPropagator::PropogateError(&DiTauConstrainedFitter::ComputeMotherLorentzVectorPar,exppar,expcov);
  return LorentzVectorParticle(m,mcov,-1.0*fabs(PDGInfo::tau_minus)*c,c,b);
}

TVectorD DiTauConstrainedFitter::Value(TVectorD &v){
  TLorentzVector Taua1,Taumu;
  double ZMass;
  CovertParToObjects(v,Taua1,Taumu,ZMass);

//     std::cout<<"DiTauConstrainedFitter  A1"<<Taua1.Px()<<"  "<<Taua1.Py()<<"  "<<Taua1.Pz()<<"  " <<Taua1.M()<<std::endl;
//     std::cout<<"DiTauConstrainedFitter  Mu"<<Taumu.Px()<<"  "<<Taumu.Py()<<"  "<<Taumu.Pz()<<"  " <<Taumu.M()<<std::endl;

  TLorentzVector z=Taua1+Taumu;
  TVectorD d(3);
//   std::cout<<"|||||||||||||||||||||||||||||| " <<std::endl;
//   std::cout<<"|||||||||||||||||||||||||||||| " <<std::endl;
//   std::cout<<"|||||||||||||||||||||||||||||| " <<std::endl;
//   std::cout<<"|||||||||||||||||||||||||||||| " <<std::endl;
//   std::cout<<"|||||||||||||||||||||||||||||| " <<std::endl;
//std::cout<<"\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ " <<std::endl;



  d(0) = z.M() - ZMass;
//  d(0) = z.M() - 91.5;
  d(1) = Taua1.Pt() - Taumu.Pt();
  d(2) = Taua1.Py() + Taumu.Py();

   std::cout<<"d(0) " << d(0)<<"  random Mass " <<ZMass <<"  ZMass   " << z.M() <<std::endl;
   std::cout<<"d(1) " << d(1)<<std::endl;
      std::cout<<"d(2) " << d(2)<<std::endl;
      std::cout<<"CSumm " << d(0) + d(1) + d(2)<<std::endl;
   //std::cout<<"CSumm " << d(0) + d(1) <<std::endl;

  return d;
}


void DiTauConstrainedFitter::CovertParToObjects(TVectorD &v,TLorentzVector &Taua1,TLorentzVector &Taumu, double &Zmass){
  Taua1=TLorentzVector(v(taua1_px),v(taua1_py),v(taua1_pz),sqrt(1.777*1.777+v(taua1_px)*v(taua1_px)+v(taua1_py)*v(taua1_py)+v(taua1_pz)*v(taua1_pz)));
  Taumu=TLorentzVector(v(taumu_px),v(taumu_py),v(taumu_pz),sqrt(1.777*1.777+v(taumu_px)*v(taumu_px)+v(taumu_py)*v(taumu_py)+v(taumu_pz)*v(taumu_pz)));
  Zmass = 91.5;
}

   
 bool DiTauConstrainedFitter::Fit(){

   std::cout<<" DiTauConstrainedFitter::Fit  "<<   std::endl;
   return LagrangeMultipliersFitter::Fit();
 }
  


LorentzVectorParticle  
DiTauConstrainedFitter::EstimateTauMu(TVector3 PV,TMatrixTSym<double>  PVCov, TVector3 SV, TMatrixTSym<double>  SVCov,TrackParticle MuTrack, LorentzVectorParticle TauA1, double ZMassR){


   TVector3 Point;
   // Define Tau direction
   TVector3 TauDir =  SV + PV;
   TVector3 TauDirError(sqrt(SVCov(0,0) + PVCov(0,0)),
			sqrt(SVCov(1,1) + PVCov(1,1)),
			sqrt(SVCov(2,2) + PVCov(2,2)));

   // Vector parameters
   double dxy   =MuTrack.Parameter(TrackParticle::dxy);
   double kappa =MuTrack.Parameter(TrackParticle::kappa);
   double phi0  =MuTrack.Parameter(TrackParticle::phi);
   double lam   =MuTrack.Parameter(TrackParticle::lambda);
   double dz    =MuTrack.Parameter(TrackParticle::dz);
   double c     =MuTrack.Charge();

   TVector3 ReferencePoint(c*dxy*sin(phi0) +PV.X() ,-c*dxy*cos(phi0)-PV.Y(),dz+PV.Z());
   TVector3 A1SV = -SV + PV;
   //------------------------------   
   //all WRT PV

   //=====
   double phiAnot  = atan2(A1SV.Y(), A1SV.X());
   double xpoca2Anot = dxy*sin(phi0) - PV.X();
   double ypoca2Anot = -dxy*cos(phi0) - PV.Y();
   double aAnot = tan(phi0);
   double bAnot = ypoca2Anot - aAnot*xpoca2Anot;
   double r = sqrt( pow(bAnot/(tan(phiAnot) - aAnot )  ,2) + pow(bAnot*tan(phiAnot)/(tan(phiAnot) - aAnot) ,2));//(bAnot)/(tan(phiAnot) - tan(phi0))/cos(phiAnot);


   double bz0 = fabs(dxy) - dz/lam ;
   double ZNeu = lam*(r +  dxy)  - dz;
   double XNeu = r*cos(phiAnot);
   double YNeu = r*sin(phiAnot);

   // std::cout<<"scalar product PV coord  "<<(xpoca2Anot +PV.X())*(xpoca2Anot - XNeu) + (ypoca2Anot+PV.Y())*(ypoca2Anot - YNeu) <<"  and phis  delta " <<fabs(atan2(YNeu ,XNeu )  - atan2(A1SV.Y(), A1SV.X()) )<<std::endl;
   std::cout<<"X,Y,Z  "<< XNeu <<"  " <<YNeu<<"  " <<ZNeu<<std::endl;
   std::cout<<"X,Y,Z  "<< XNeu <<"  " <<YNeu<<"  " <<lam*(r +  dxy)  + dz<<std::endl;
 
   //-------------------------------
   //------------------------------   
   //all WRT 000

   double phiAnot1  = atan2(A1SV.Y(), A1SV.X());

    double xpoca2Anot1 = dxy*sin(phi0);
    double ypoca2Anot1 = -dxy*cos(phi0);
    double at = (SV-PV).Y()/(SV-PV).X(), bt =SV.Y() - at*SV.X() ;
    double am = tan(phi0), bm =ypoca2Anot1-am*xpoca2Anot1;

    double XNeu1 = (bt-bm)/(am-at);
    double YNeu1 = am*XNeu1 + bm;

    double r1 = sqrt(XNeu1*XNeu1  + YNeu1*YNeu1  );
    double phiWRT000 = atan2(YNeu1,XNeu1);

    double bz01 = fabs(dxy) - dz/lam ;
    double ZNeu1 = lam*(r1 +  dxy)  - dz;

//     std::cout<<"scalar product 000 coord  "<<xpoca2Anot1*(xpoca2Anot1 - XNeu1) + ypoca2Anot1*(ypoca2Anot1 - YNeu1) <<"  and phis  delta " <<fabs(atan2((YNeu1 - PV.Y()),(XNeu1 - PV.X()) )  - atan2(A1SV.Y(), A1SV.X()) )<<std::endl;
//    std::cout<<"check coord "<<XNeu1<<"  " <<YNeu1<<" =   "<< r1*cos(phiWRT000)<< "  "<<r1*sin(phiWRT000) <<std::endl;
//    std::cout<<"transformation  "<< XNeu+PV.X()<< " ==  "<<XNeu1<<"  "<< YNeu+ PV.Y()<<"  == " <<YNeu1 <<std::endl;





//    //=====

//    //-------------------------------



//    //    double ZNeu = lam*(r +  dxy)  - dz; // + PV.Z();
//    //    double XNeu = r*cos(phiAnot); //+ PV.X();
//    //    double YNeu = r*sin(phiAnot);// + PV.Y();
//    std::cout<<"XNeu1  "<<  XNeu1 <<" YNeu  "<< YNeu1<<" ZNeu1  "<< ZNeu1<<std::endl;
//    std::cout<<"xpoca2Anot1   "<<  xpoca2Anot1 <<" ypoca2Anot   "<<ypoca2Anot1<<std::endl;



//    std::cout<<"XNeu  "<<  XNeu  <<" YNeu  "<< YNeu<<" ZNeu  "<< ZNeu<<std::endl;
//    std::cout<<"xpoca2Anot  "<<  xpoca2Anot <<" ypoca2Anot   "<<ypoca2Anot <<std::endl;


//    std::cout<<"PVx  "<< PV.X()   <<" PVy  "<<PV.Y()   <<" PVz  "<< PV.Z()  <<std::endl;
//    std::cout<<"SVx  "<< SV.X()   <<" SVy  "<<SV.Y()   <<" SVz  "<< SV.Z()  <<std::endl;
 
   //    std::cout<<"TauDirError  "<<  TauDirError.X() <<"   "<< TauDirError.Y() << "   " << TauDirError.Z()<<std::endl;
   //    std::cout<<"TauDirr  = "<<  A1SV.Y() <<" /  "<< A1SV.X() << "  =   " << A1SV.Y()/A1SV.X()<< " deltaTanPhi " <<sqrt(TauDirError.Y()*TauDirError.Y()/A1SV.X()/A1SV.X()   + TauDirError.X()*TauDirError.X()*A1SV.Y()*A1SV.Y()/A1SV.X()/A1SV.X()/A1SV.X()/A1SV.X() )  <<"  deltaPhi " << sqrt(TauDirError.Y()*TauDirError.Y()/A1SV.X()/A1SV.X()   + TauDirError.X()*TauDirError.X()*A1SV.Y()*A1SV.Y()/A1SV.X()/A1SV.X()/A1SV.X()/A1SV.X() ) *cos(phiAnot)*cos(phiAnot)<<std::endl;



   //---- Covariance of original parameters (dxy, phi0,dz,lambda)
   TMatrixT<double>  HelixCov;
   HelixCov.ResizeTo(5,5);
   
   HelixCov(0,0) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy);
   HelixCov(0,1) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::phi);
   HelixCov(0,2) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dz);
   HelixCov(0,3) = MuTrack.Covariance(TrackParticle::dxy,TrackParticle::lambda);
   HelixCov(0,4) = 0;

   HelixCov(1,0) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dxy);
   HelixCov(1,1) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi);
   HelixCov(1,2) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::dz);
   HelixCov(1,3) = MuTrack.Covariance(TrackParticle::phi,TrackParticle::lambda);
   HelixCov(1,4) = 0;

   HelixCov(2,0) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dxy);
   HelixCov(2,1) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::phi);
   HelixCov(2,2) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz);
   HelixCov(2,3) = MuTrack.Covariance(TrackParticle::dz,TrackParticle::lambda);
   HelixCov(2,4) = 0;

   HelixCov(3,0) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dxy);
   HelixCov(3,1) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::phi);
   HelixCov(3,2) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::dz);
   HelixCov(3,3) = MuTrack.Covariance(TrackParticle::lambda,TrackParticle::lambda);
   HelixCov(3,4) = 0;

   HelixCov(4,0) = 0;
   HelixCov(4,1) = 0;
   HelixCov(4,2) = 0;
   HelixCov(4,3) = 0;
   HelixCov(4,4) = sqrt(TauDirError.Y()*TauDirError.Y()/A1SV.X()/A1SV.X()   + TauDirError.X()*TauDirError.X()*A1SV.Y()*A1SV.Y()/A1SV.X()/A1SV.X()/A1SV.X()/A1SV.X() ) *cos(phiAnot)*cos(phiAnot);
   
   //   std::cout<<"DiTauConstraint  HelixCov"<<std::endl;
   //   for(int str =0; str < HelixCov.GetNrows(); str++){
   //     for(int kol =0; kol < HelixCov.GetNcols(); kol++){
   //       std::cout<<"  "<< HelixCov(str,kol)<<"  ";
   //     }    
   //     std::cout<<std::endl;
   //   }

  //-------------------
  double drdd = 2*sin(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  double drdphi0 = 2*dxy*cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))   - 2*dxy*sin(phi0)*cos(phiAnot)/cos(phi0)/cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  double drdphi =2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  TMatrixT<double>  DerivativesXYZToHelix;

  DerivativesXYZToHelix.ResizeTo(3,5);




  // Set derivatives matrix  (dxy,phi0,lambda,z,phi)
  //  std::cout<<" PARAMETERES:    drdd  "<< drdd <<" drdphi0  "<<drdphi0<<"  drdphi  "<<drdphi<< " lam  "<< lam <<std::endl;
  DerivativesXYZToHelix(0,0)  = drdd*cos(phiAnot);
  DerivativesXYZToHelix(0,1)  = drdphi0*cos(phiAnot);
  DerivativesXYZToHelix(0,2)  = 0;
  DerivativesXYZToHelix(0,3)  = 0;
  DerivativesXYZToHelix(0,4)  = drdphi*cos(phiAnot) - r*sin(phiAnot);

  DerivativesXYZToHelix(1,0)  = drdd*sin(phiAnot);
  DerivativesXYZToHelix(1,1)  = drdphi0*sin(phiAnot);
  DerivativesXYZToHelix(1,2)  = 0;
  DerivativesXYZToHelix(1,3)  = 0;
  DerivativesXYZToHelix(1,4)  = drdphi*sin(phiAnot) - r*cos(phiAnot);

  DerivativesXYZToHelix(2,0)  = lam;
  DerivativesXYZToHelix(2,1)  = lam*drdphi0;
  DerivativesXYZToHelix(2,2)  = r+dxy;
  DerivativesXYZToHelix(2,3)  = -1;
  DerivativesXYZToHelix(2,4)  = lam*drdphi;

  TVector3 NormalisedTauMuDirection;  TVector3 NormalisedTauMuDirectionWRT000;
  TVector3 NormalisedTauA1Direction;


  //  std::cout<<"transformation  "<< XNeu+PV.X()<< " ==  "<<XNeu1<<"  "<< YNeu+ PV.Y()<<"  == " <<YNeu1 <<std::endl;

  TVector3 PointWRTPV(XNeu,YNeu,ZNeu);
  TVector3 PointWRT000(XNeu1 - PV.X(),YNeu1 - PV.Y(),ZNeu1 - PV.Z());



  NormalisedTauMuDirection.SetX(PointWRTPV.X()/PointWRTPV.Mag());
  NormalisedTauMuDirection.SetY(PointWRTPV.Y()/PointWRTPV.Mag());
  NormalisedTauMuDirection.SetZ(PointWRTPV.Z()/PointWRTPV.Mag());
  // wrt to 000
  //    NormalisedTauMuDirection.SetX(PointWRT000.X()/PointWRTPV.Mag());
  //    NormalisedTauMuDirection.SetY(PointWRT000.Y()/PointWRTPV.Mag());
  //    NormalisedTauMuDirection.SetZ(PointWRT000.Z()/PointWRTPV.Mag());



  NormalisedTauA1Direction.SetX(TauA1.LV().Px()/TauA1.LV().P());
  NormalisedTauA1Direction.SetY(TauA1.LV().Py()/TauA1.LV().P());
  NormalisedTauA1Direction.SetZ(TauA1.LV().Pz()/TauA1.LV().P());
  

  //  std::cout<<"chek new phi  "<<  atan2(YNeu,XNeu) <<" SV  "<< phiAnot<<" fabs(new Phi)   " <<fabs(phiAnot - atan2(YNeu,XNeu))<<"  " <<atan2(NormalisedTauMuDirection.Y(),  NormalisedTauMuDirection.X() )<<std::endl;

  double cosTauTau = NormalisedTauMuDirection.X()* NormalisedTauA1Direction.X() + NormalisedTauMuDirection.Y()* NormalisedTauA1Direction.Y() + NormalisedTauMuDirection.Z()* NormalisedTauA1Direction.Z();


  TMatrixT<double> DerivativesXYZToHelixT=DerivativesXYZToHelix; DerivativesXYZToHelixT.T();
  TMatrixT<double> CovXYZFrame=DerivativesXYZToHelix*HelixCov*DerivativesXYZToHelixT;
  
  //  double TauA1deltaP = sqrt(   (pow(TauA1.LV().Px()*TauA1.Covariance(3,3) ,2)  +  pow(TauA1.LV().Py()*TauA1.Covariance(4,4),2)  +  pow(TauA1.LV().Pz()*TauA1.Covariance(5,5),2)    )/TauA1.LV().P()/TauA1.LV().P()  )  ;
  double TauA1deltaP = sqrt(   (pow(TauA1.LV().Px(),2)*TauA1.Covariance(3,3)   +  pow(TauA1.LV().Py(),2)*TauA1.Covariance(4,4)  +  pow(TauA1.LV().Pz(),2)*TauA1.Covariance(5,5)    )/TauA1.LV().P()/TauA1.LV().P()  )  ;
  double TauMudeltaP = TauA1deltaP*pow(ZMassR/2/TauA1.LV().P(),2)/(1-cosTauTau);
  double TauMuP = ZMassR*ZMassR/2/(1-cosTauTau)/TauA1.LV().P();
                 
  TLorentzVector TauMuEstimate;
  TauMuEstimate.SetXYZM(TauMuP*NormalisedTauMuDirection.X(),TauMuP*NormalisedTauMuDirection.Y(),TauMuP*NormalisedTauMuDirection.Z(),1.777);

  std::cout<<" mass check  "<<ZMassR <<" INV  "<< (TauMuEstimate+TauA1.LV()).M()<<" TauMuP   " << TauMuP<<std::endl;

  std::cout<<"   costTau   "<<cosTauTau <<" mass equation   " << 2*TauA1.LV().P()*(1-cosTauTau)* TauMuP<<std::endl ;

  
  std::cout<<"TauA1 covariance "<< sqrt(TauA1.Covariance(3,3))<<" "<<sqrt(TauA1.Covariance(4,4))<<" "<<sqrt(TauA1.Covariance(5,5))<<" TauA1 deltaPa  " <<TauA1deltaP <<" TauMu deltaPa  " <<TauMudeltaP<<std::endl;
//   std::cout<<" TauMuEstimate  "<< TauMuEstimate.Px()<<" "<<TauMuEstimate.Py()<<"  "<<TauMuEstimate.Pz()<<std::endl;
//   std::cout<<" TauA1P  "<< TauA1.LV().Px()<<" "<<TauA1.LV().Py()<<"  "<<TauA1.LV().Pz()<<std::endl;
//   std::cout<<" deltaPhi Mu - TauA1  "<< fabs(TauA1.LV().Phi() - TauMuEstimate.Phi())<<std::endl;
//   std::cout<<" deltaPhi TauMu - TauA1  "<< TauA1.LV().Phi() - phi0  <<std::endl;
//   std::cout<<" DiTauMass  "<< (TauA1.LV() + TauMuEstimate).M()  <<" generated value  " <<ZMassR <<std::endl;



   //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   Compute Covariance in terms of angles only a

  double br = tan(phi0)*cos(phiAnot) - sin(phiAnot);
  double Zc = ZNeu;

  double MinusSintheta = -sqrt(1  - Zc*Zc/(r*r + Zc*Zc) );

  double cosTheta = Zc/sqrt(r*r + Zc*Zc) ;
  double sinTheta = -MinusSintheta;
  if(fabs(phiAnot - TauA1.LV().Phi()) < 0.05) phiAnot = phiAnot - TMath::Pi();
  std::cout<<" deltaPhi's " <<fabs(phiAnot - TauA1.LV().Phi())<<std::endl;

  std::cout<<"   costTauTau2   "<< cos(acos(cosTheta) + TauA1.LV().Theta())  <<std::endl;
  //-derivaticves

  double drzdd = 2*sin(phi0)/br/Zc  - 2*dxy*sin(phi0)*lam/br/Zc/Zc;
  double drzdphi0 = (2*dxy*cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))   - 2*dxy*sin(phi0)*cos(phiAnot)/cos(phi0)/cos(phi0)/(tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot)))/Zc;

  double drzdlam = -2*dxy*sin(phi0)*(r+dxy)/br/Zc/Zc;
  double drzdz0 = -2*dxy*sin(phi0)/br/Zc;

  double drzdphi = drdphi*(1 - r*lam/Zc)/Zc;//2*dxy*sin(phi0)*(cos(phiAnot) + tan(phi0)*sin(phiAnot))/ (tan(phi0)*cos(phiAnot) - sin(phiAnot))/(tan(phi0)*cos(phiAnot) - sin(phiAnot));
  double dcosThetadrz = -r/Zc/sqrt(pow(1 + r*r/Zc/Zc,3));

  double dThetadPhi = -dcosThetadrz*drzdphi/sinTheta;
  //-derivaticves

  double cosTauTau2 = cos(TauA1.LV().Theta() + acos(cosTheta));
  double sinTauTau2 = sin(TauA1.LV().Theta() + acos(cosTheta));


  double TauMudeltaP_2 = TauA1deltaP*pow(ZMassR/2/TauA1.LV().P(),2)/(1-cosTauTau2);

  double TauMuP_2 =  ZMassR*ZMassR/2/(1-cosTauTau2)/TauA1.LV().P();

  double dPdTheta  = ZMassR*ZMassR/2/TauA1.LV().P()*sinTauTau2/(1-cosTauTau2)/(1-cosTauTau2);
  double dPdPhi  = dPdTheta*dThetadPhi;


  std::cout<<"TauMuP2  pm TauMudeltaP2 -----  "<< TauMuP_2<<"  " << TauMudeltaP_2<<" TauMuP_2  " <<TauMuP_2 <<std::endl;
 
  std::cout<<"TauMuP2  pm TauMudeltaP2 -----  "<< TauMuP_2<<"  " << TauMudeltaP_2<<" TauMuP_2  " <<TauMuP_2 <<std::endl;


  TLorentzVector TauMuEstimate2;
  TauMuEstimate2.SetXYZM(TauMuP_2*cos(phiAnot)*sinTheta,TauMuP_2*sin(phiAnot)*sinTheta,TauMuP_2*cosTheta,1.777);


  std::cout<<"check kinematic 1 TauMuP  "<< TauMuP<<"   1- cos " << 1-cosTauTau<<"   TauA1P  " <<TauA1.LV().P() << " Mass   " << (TauA1.LV() +TauMuEstimate).M()   <<std::endl;
  std::cout<<"check kinematic 2 TauMuP  "<< TauMuP_2 <<"  1- cs" <<1-cosTauTau2 <<"  TauA1P  " <<TauA1.LV().P() <<"  Mass " << (TauA1.LV() +TauMuEstimate2).M() <<std::endl;


  std::cout<<"compare phis   "<< TauMuEstimate2.Phi()  <<" " << TauA1.LV().Phi()<<std::endl;
 
  std::cout<<"TauMuP_2*TauA1.LV().P()*(1-cosTauTau2)  "<< TauMuP_2*TauA1.LV().P()*cosTauTau2 <<" TauMuEstimate2.E " << TauMuEstimate2.Px()*TauA1.LV().Px() + TauMuEstimate2.Py()*TauA1.LV().Py()+ TauMuEstimate2.Pz()*TauA1.LV().Pz() << " check 3 " << TauMuP_2*TauA1.LV().P()*cos(TauMuEstimate2.Theta() +TauA1.LV().Theta() )<<std::endl;





   std::cout<<"LorentzVectorParticle 1 "<< TauA1.LV().Px()  <<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<< TauA1.LV().Py()   <<std::endl;
   std::cout<<"LorentzVectorParticle 1 "<< TauA1.LV().Pz()   <<std::endl;
   
   std::cout<<"LorentzVectorParticle 2 "<< TauMuEstimate2.Px()   <<std::endl;
   std::cout<<"LorentzVectorParticle 2 "<<TauMuEstimate2.Py()   <<std::endl;
   std::cout<<"LorentzVectorParticle 2 "<< TauMuEstimate2.Pz()   <<std::endl;
   std::cout<<" mass check  "<<ZMassR <<" INV  "<< (TauMuEstimate2+TauA1.LV()).M()<<std::endl;



   std::cout<<"check delta  1 "<< (TauA1.LV()   +  TauMuEstimate2).M()<<std::endl;
   std::cout<<"check delta  1 "<< TauA1.LV().Px()  + TauMuEstimate2.Px() <<std::endl;
   std::cout<<"check delta  1 "<< TauA1.LV().Py()  + TauMuEstimate2.Py()  <<std::endl;



  //double dPd

  TMatrixT<double>  DerivativesHelixToAngles;
  DerivativesHelixToAngles.ResizeTo(2,5);




  // Set derivatives matrix  (dxy,phi0,lambda,z,phi)
  //  std::cout<<" PARAMETERES:    drdd  "<< drdd <<" drdphi0  "<<drdphi0<<"  drdphi  "<<drdphi<< " lam  "<< lam <<std::endl;
  DerivativesHelixToAngles(0,0)  = 0;
  DerivativesHelixToAngles(0,1)  = 0;
  DerivativesHelixToAngles(0,2)  = 0;
  DerivativesHelixToAngles(0,3)  = 0;
  DerivativesHelixToAngles(0,4)  = 1;

  DerivativesHelixToAngles(1,0)  = dcosThetadrz*drzdd/MinusSintheta;
  DerivativesHelixToAngles(1,1)  = dcosThetadrz*drzdphi0/MinusSintheta;
  DerivativesHelixToAngles(1,2)  = dcosThetadrz*drzdlam/MinusSintheta;
  DerivativesHelixToAngles(1,3)  = dcosThetadrz*drzdz0/MinusSintheta;
  DerivativesHelixToAngles(1,4)  = dcosThetadrz*drzdphi/MinusSintheta;

  TMatrixT<double> DerivativesHelixToAnglesT=DerivativesHelixToAngles; DerivativesHelixToAnglesT.T();
  TMatrixT<double> CovAngleFrame=DerivativesHelixToAngles*HelixCov*DerivativesHelixToAnglesT;

  TMatrixT<double> CovPPhiThetaFrame;
  CovPPhiThetaFrame.ResizeTo(3,3);
//---------- dp dphi dtheta
      CovPPhiThetaFrame(0,0) = TauMudeltaP;
      CovPPhiThetaFrame(0,1) = 0;
      CovPPhiThetaFrame(0,2) = 0;

      CovPPhiThetaFrame(1,0) = 0;
      CovPPhiThetaFrame(1,1) = CovAngleFrame(0,0);
      CovPPhiThetaFrame(1,2) = CovAngleFrame(0,1);

      CovPPhiThetaFrame(2,0) = 0;
      CovPPhiThetaFrame(2,1) = CovAngleFrame(1,0);
      CovPPhiThetaFrame(2,2) = CovAngleFrame(1,1);
//    include correlations

//      CovPPhiThetaFrame(0,0) = TauMudeltaP;
//      CovPPhiThetaFrame(0,1) = dPdPhi*dPdPhi*CovAngleFrame(0,0);
//      CovPPhiThetaFrame(0,2) = dPdTheta*dPdTheta*CovAngleFrame(1,1);

//      CovPPhiThetaFrame(1,0) = dPdPhi*dPdPhi*CovAngleFrame(0,0);
//      CovPPhiThetaFrame(1,1) = CovAngleFrame(0,0);
//      CovPPhiThetaFrame(1,2) = CovAngleFrame(0,1);

//      CovPPhiThetaFrame(2,0) = dPdTheta*dPdTheta*CovAngleFrame(1,1);
//      CovPPhiThetaFrame(2,1) = CovAngleFrame(1,0);
//      CovPPhiThetaFrame(2,2) = CovAngleFrame(1,1);


   std::cout<<" derivatives   "<< " dPdTheta  "<<dPdTheta*dPdTheta*CovAngleFrame(1,1)<<"  dPdPhi  "<<dPdPhi*dPdPhi*CovAngleFrame(0,0)<<std::endl;



   std::cout<<"CovAngleFrame -----"<<std::endl;
   for(int str =0; str < CovAngleFrame.GetNrows(); str++){
     for(int kol =0; kol < CovAngleFrame.GetNcols(); kol++){
       std::cout<<"  "<< CovAngleFrame(str,kol)<<"  ";
     }    
     std::cout<<std::endl;
   }
   std::cout<<"CovPPhiThetaFrame -----"<<std::endl;
   for(int str =0; str < CovPPhiThetaFrame.GetNrows(); str++){
     for(int kol =0; kol < CovPPhiThetaFrame.GetNcols(); kol++){
       std::cout<<"  "<< CovPPhiThetaFrame(str,kol)<<"  ";
     }    
     std::cout<<std::endl;
   }


  TMatrixT<double> DirevativesPPhiThetaToPxPyPz;
  DirevativesPPhiThetaToPxPyPz.ResizeTo(3,3);

  DirevativesPPhiThetaToPxPyPz(0,0) =  cos(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(0,1) = -TauMuP_2*sin(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(0,2) =  TauMuP_2*cos(phiAnot)*cosTheta;


  DirevativesPPhiThetaToPxPyPz(1,0) = sin(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(1,1) = TauMuP_2*cos(phiAnot)*sinTheta;
  DirevativesPPhiThetaToPxPyPz(1,2) = TauMuP_2*sin(phiAnot)*cosTheta;

  DirevativesPPhiThetaToPxPyPz(2,0) = cosTheta;
  DirevativesPPhiThetaToPxPyPz(2,1) = 0;
  DirevativesPPhiThetaToPxPyPz(2,2) =-TauMuP_2*sinTheta;



  TMatrixT<double> DirevativesPPhiThetaToPxPyPzT=DirevativesPPhiThetaToPxPyPz; DirevativesPPhiThetaToPxPyPzT.T();
  TMatrixT<double> CovPxPyPzFrame=DirevativesPPhiThetaToPxPyPz*CovPPhiThetaFrame*DirevativesPPhiThetaToPxPyPzT;


  std::cout<<"TauMuP  pm TauMudeltaP -----  "<< TauMuP<<"  " << TauMudeltaP<<std::endl;

  std::cout<<"DirevativesPPhiThetaToPxPyPz -----"<<std::endl;
  for(int str =0; str < DirevativesPPhiThetaToPxPyPz.GetNrows(); str++){
    for(int kol =0; kol < DirevativesPPhiThetaToPxPyPz.GetNcols(); kol++){
      std::cout<<"  "<< DirevativesPPhiThetaToPxPyPz(str,kol)<<"  ";
    }    
    std::cout<<std::endl;
  }


  std::cout<<"CovPxPyPzFrame -----"<<std::endl;
  for(int str =0; str < CovPxPyPzFrame.GetNrows(); str++){
    for(int kol =0; kol < CovPxPyPzFrame.GetNcols(); kol++){
      std::cout<<"  "<< CovPxPyPzFrame(str,kol)<<"  ";
    }    
    std::cout<<std::endl;
  }


   //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Compute Covariance in terms of angles only a


  TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
  TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);

  par(LorentzVectorParticle::vx,0)=XNeu;
  par(LorentzVectorParticle::vy,0)=YNeu;
  par(LorentzVectorParticle::vz,0)=ZNeu;
  par(LorentzVectorParticle::px,0)=TauMuEstimate2.Px();
  par(LorentzVectorParticle::py,0)=TauMuEstimate2.Py();
  par(LorentzVectorParticle::pz,0)=TauMuEstimate2.Pz();
  par(LorentzVectorParticle::m,0) =1.777;


  //---- FillVertexCov
  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      //     Cov(i+3,j+3)=TauA1deltaP + Cov(i+3,j+3);//TauMuEstimate.P()*CovXYZFrame(i,j);
      Cov(i,j)=CovXYZFrame(i,j);

      //    Cov(i,j)=CovPxPyPzFrame(i,j);

      //      Cov(i+3,j+3)=TauA1.Covariance(3,3) + Cov(i+3,j+3);//TauMuEstimate.P()*CovXYZFrame(i,j);
      
    }
  }
  //---- FillVertexCov


  TMatrixT<double>    VectorEForCovariance;
  TMatrixT<double>    VectorRForCovariance;

  TMatrixT<double>    MatrixEForCovariance;
  TMatrixT<double>    FinCovariance;


  VectorEForCovariance.ResizeTo(3,1);
  VectorRForCovariance.ResizeTo(3,1);

  MatrixEForCovariance.ResizeTo(3,3);
  FinCovariance.ResizeTo(3,3);

  VectorEForCovariance(0,0) = TauMuEstimate.P();
  VectorEForCovariance(0,1) = TauMuEstimate.P();
  VectorEForCovariance(0,2) = TauMuEstimate.P();

  VectorRForCovariance(0,0) = XNeu;
  VectorRForCovariance(0,1) = YNeu;
  VectorRForCovariance(0,2) = ZNeu;

  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){
      MatrixEForCovariance(i,j)=0;
      if(i==j)MatrixEForCovariance(i,j) = TauMudeltaP*(pow(XNeu,2)+pow(YNeu,2)+pow(ZNeu,2) ); 
    }
  }
  FinCovariance = MatrixEForCovariance + TauMuEstimate.P()*TauMuEstimate.P()*CovXYZFrame;


  for(int i=0; i<LorentzVectorParticle::NVertex; i++){
    for(int j=0; j<LorentzVectorParticle::NVertex; j++){


      //        Cov(i+3,j+3)=TauMuEstimate.P()*CovXYZFrame(i,j);

      // Cov(i+3,j+3)=FinCovariance(i,j);

      Cov(i+3,j+3)=CovPxPyPzFrame(i,j);
      
    }
  }//Cov(i,j)=CovPxPyPzFrame(i,j);

        std::cout<<"FinCovariance "<<std::endl;
        for(int str =0; str < FinCovariance.GetNrows(); str++){
         for(int kol =0; kol <FinCovariance.GetNcols(); kol++){
            std::cout<<"  "<< FinCovariance(str,kol)<<"  ";
          }    
        std::cout<<std::endl;
        }


        std::cout<<"PointWRTPV X-----"<<PointWRTPV.X()<<"  "<<std::endl;
        std::cout<<"PointWRTPV Y-----"<<PointWRTPV.Y()<<"  "<<std::endl;
        std::cout<<"PointWRTPV Z-----"<<PointWRTPV.Z()<<"  "<<std::endl;


        std::cout<<"DerivativesXYZToHelix -----"<<std::endl;
        for(int str =0; str < DerivativesXYZToHelix.GetNrows(); str++){
          for(int kol =0; kol < DerivativesXYZToHelix.GetNcols(); kol++){
            std::cout<<"  "<< DerivativesXYZToHelix(str,kol)<<"  ";
          }    
        std::cout<<std::endl;
        }





        std::cout<<"HelixCove -----"<<std::endl;
        for(int str =0; str < HelixCov.GetNrows(); str++){
          for(int kol =0; kol < HelixCov.GetNcols(); kol++){
            std::cout<<"  "<< HelixCov(str,kol)<<"  ";
          }    
        std::cout<<std::endl;
        }





       std::cout<<"DerivativesXYZToHelix  transpose  "<<std::endl;
       for(int str =0; str < DerivativesXYZToHelixT.GetNrows(); str++){
        for(int kol =0; kol < DerivativesXYZToHelixT.GetNcols(); kol++){
           std::cout<<"  "<< DerivativesXYZToHelixT(str,kol)<<"  ";
         }    
       std::cout<<std::endl;
       }




     std::cout<<"CovXYZFrame "<<std::endl;
     for(int str =0; str < CovXYZFrame.GetNrows(); str++){
       for(int kol =0; kol < CovXYZFrame.GetNcols(); kol++){
         std::cout<<"     "<< CovXYZFrame(str,kol)<<"     ";
       }    
     std::cout<<std::endl;
     }

    return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);


 //   //------------------------------   

   //    TVector3 ReferencePoint(1,2,0);
   //    TVector3 A1SV(-2,-1,0);


//    double xpoca2 = dxy*sin(phi0);
//    double ypoca2 = -dxy*cos(phi0);

//    double aMuonTrack  = tan(MuonLV.Phi());
//    double bMuonTrack = ReferencePoint.Y() - aMuonTrack*ReferencePoint.X();


//    double aTau = (SV.Y() - PV.Y())/(SV.X() - PV.X());
//    double bTau = (SV.X()*PV.Y() - PV.X()*SV.Y())/(SV.X() - PV.X());


//    double xTest =  0.7;
//    double yTest = aTau*xTest + bTau;


//    TVector3 A1TauDirection(SV.X() - PV.X(),SV.Y() - PV.Y(),SV.Z() - PV.Z()); // With respect to PV

//  //   std::cout<<"opposite to tau "<< tan(A1TauDirection.Phi() - 3.1415)<< "  "<<tan(A1TauDirection.Phi() + 3.1415)<<std::endl;

//    double a = tan(A1TauDirection.Phi()); 
// //    double InterX =(bMuonTrack-bTau)/(aTau-aMuonTrack);
// //    double InterY = (bMuonTrack*aTau-aMuonTrack*bTau)/(aTau-aMuonTrack);

//    double phi1 =atan2(A1SV.Y(), A1SV.X());
//    double phi =phi1;
// //    if(phi1 > 0 && phi1 <  TMath::Pi())  phi  = phi1 - TMath::Pi();
// //    if(phi1 < 0 && phi1 > -TMath::Pi())  phi  = phi1 + TMath::Pi();

//    double InterX = (bMuonTrack - bTau)/(aTau - aMuonTrack);//( PV.Y() - bMuonTrack - PV.X()*tan(phi))/(aMuonTrack  - tan(phi));
//    double InterY = aMuonTrack*InterX + bMuonTrack;
//    double RotatedPhi =0;








//    if(fabs(atan2(SV.Y()- InterY,SV.X()- InterX ) - atan2(SV.Y()-PV.Y(),SV.X()-PV.X())) > 0.1 and atan2(SV.Y()-PV.Y(),SV.X()-PV.X()) < 0 )  RotatedPhi = phi + TMath::Pi();
//    if(fabs(atan2(SV.Y()- InterY,SV.X()- InterX ) - atan2(SV.Y()-PV.Y(),SV.X()-PV.X())) > 0.1 and atan2(SV.Y()-PV.Y(),SV.X()-PV.X()) > 0 )  RotatedPhi = phi - TMath::Pi();

//    double InterXFin = ( PV.Y() - bMuonTrack - PV.X()*tan(RotatedPhi))/(aMuonTrack  - tan(RotatedPhi));
//    double InterYFin = aMuonTrack*InterX + bMuonTrack;


//    TVector3  DirectionWRTPV( InterX,InterY ,0);   
//    TVector3  Direction2WRTPV(A1SV.X() ,A1SV.Y() ,0);

//    std::cout<<"PV  "<<  PV.X() <<" "<< PV.Y()<<std::endl;
//    std::cout<<"SV  "<<  SV.X() <<" "<< SV.Y()<<std::endl;
//    std::cout<<"POCA  "<<  ReferencePoint.X() <<" "<< ReferencePoint.Y()<<std::endl;

//    std::cout<<"Point  "<<  InterX <<" "<< InterY<<std::endl;
//    std::cout<<"PointFin  "<<  InterXFin <<" "<<InterYFin <<std::endl;
//    std::cout<<"Test  "<<  xTest <<" "<< yTest<<std::endl;
//    std::cout<<"MuonPhi  "<<  MuonLV.Phi() <<" "<<  atan2(-InterY+0.001,-InterX+0.001)<<"  "<<phi0 <<std::endl;


//    std::cout<<"comparephis  "<<  atan2(SV.Y()-PV.Y(),SV.X()-PV.X()) <<"  == "<<  atan2(SV.Y()- InterYFin,SV.X()- InterXFin )<<"  "<< atan2(SV.Y()-yTest,SV.X()- xTest ) << "fabs "<< fabs(atan2(SV.Y()-PV.Y(),SV.X()-PV.X())-atan2(SV.Y()- InterY,SV.X()- InterX ))<<std::endl;
//    std::cout<<"check that point satisfies equations  "<<  InterY << " ==   "  << aTau*InterX + bTau<<std::endl;
//    std::cout<<"check that point satisfies equations2  "<<  InterY << " ==   "  << aMuonTrack*InterX + bMuonTrack<<std::endl;




//    std::cout<<"IntersectionPointLinearApproximation DirectionWRTPV "<< DirectionWRTPV.Phi() <<" " << DirectionWRTPV.X() <<" "<< DirectionWRTPV.Y()<<std::endl;
//    std::cout<<"IntersectionPointLinearApproximation DirectionWRTPV "<< Direction2WRTPV.Phi() <<" " << Direction2WRTPV.X() <<" "<< Direction2WRTPV.Y()<<std::endl;
//    std::cout<<"POCA  "<< MuonPoca.Phi() <<" " << MuonPoca.X() <<" "<< MuonPoca.Y()<<std::endl;
 
//    Point.SetX(InterX);
//    Point.SetY(InterY);
//    Point.SetZ(0);




   
 }























///--------------------  not debugged helix approximation of track ------------




// LorentzVectorParticle DiTauConstrainedFitter::EstimateTauMu(LorentzVectorParticle Mu, double TauMuEnergyEstimate,TVector3 PVertex, TMatrixTSym<double> PVertexCov){
// //   std::cout<<"deb 19"<<std::endl; 
//   TVector3 MuonImpact = Mu.Vertex();
//   TVector3 TauMuDirEstimate  =  MuonImpact - PVertex;
//   double R = TauMuDirEstimate.Mag();

//   double phi = TauMuDirEstimate.Phi();
//   double theta =TauMuDirEstimate.Theta();

//   TMatrixTSym<double> TauMuImpactCov =  Mu.VertexCov();
//   TMatrixTSym<double> TauMuImpactPVCov;

//   TauMuImpactPVCov.ResizeTo(Mu.VertexCov().GetNrows(),Mu.VertexCov().GetNcols()); 


//   for(int i=0; i<LorentzVectorParticle::NVertex;i++){
//     for(int j=0; j<LorentzVectorParticle::NVertex;j++){
//       TauMuImpactPVCov(i,j) = TauMuEnergyEstimate*(TauMuImpactCov(i,j) + PVertexCov(i,j))/R;
//     }
//   }


//   TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
//   par(LorentzVectorParticle::vx,0)=Mu.Parameter(LorentzVectorParticle::vx);
//   par(LorentzVectorParticle::vy,0)=Mu.Parameter(LorentzVectorParticle::vy);
//   par(LorentzVectorParticle::vz,0)=Mu.Parameter(LorentzVectorParticle::vz);
//   par(LorentzVectorParticle::px,0)=TauMuEnergyEstimate*cos(phi)*sin(theta); 
//   par(LorentzVectorParticle::py,0)=TauMuEnergyEstimate*sin(phi)*sin(theta);
//   par(LorentzVectorParticle::pz,0)=TauMuEnergyEstimate*cos(theta);
//   par(LorentzVectorParticle::m,0) =1.777;

//   TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);
//   TMatrixTSym<double> pvCov=TauMuImpactPVCov;
 
//   for(int i=0; i<3; i++){
//     for(int j=0; j<3; j++){
//       if(i<LorentzVectorParticle::NVertex) Cov(i,j)=pvCov(i,j);
//       else Cov(i,j)=0;
//     }

// //     //if(i>2 && i < LorentzVectorParticle::NLorentzandVertexPar)
// //     double v=0;
// //     if(i==LorentzVectorParticle::px || i==LorentzVectorParticle::py || i==LorentzVectorParticle::pz) v=10*par(i,0)*par(i,0);
// //     if(v<1000.0) v=1000.0; // try lowing to test impact
// //     Cov(i,i)+=v;
//   }



//   for(int i=3; i<6; i++){
//     for(int j=3; j<6; j++){
//       Cov(i,j) = TauMuEnergyEstimate*TauMuEnergyEstimate*Mu.Covariance(i-3,j-3);
// //       std::cout<<" Mu.Covariance(i-3,j-3) "<< Mu.Covariance(i-3,j-3)<<"  "<<std::endl; 
// //       std::cout<<" TauMuEnergyEstimate "<< TauMuEnergyEstimate<<"  "<<std::endl; 
//     }
//   } 


//   //   std::cout<<"EstimateTauMu--->"<<std::endl; 
//   //----
 
//   for(int str =0; str < Cov.GetNrows(); str++){
//     for(int kol =0; kol < Cov.GetNcols(); kol++){
//        std::cout<<"  "<< Cov(str,kol)<<"  ";
//     }    
//      std::cout<<std::endl; 
//   }

// //   std::cout<<"deb 20"<<std::endl; 
//   return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,Mu.BField());
// }

 
// LorentzVectorParticle DiTauConstrainedFitter::ConvertTrackParticleToLorentzVectorParticle(TrackParticle MuTrack){

// //   std::cout<<"deb 211"<<std::endl; 
// //   std::cout<<"deb 212"<<std::endl; 
//   TVector3  MuonImpact;
// //   std::cout<<"deb 213"<<std::endl; 
//   double dxy = MuTrack.Parameter(TrackParticle::dxy);
// //   std::cout<<"deb 214"<<std::endl; 
//   double dz = MuTrack.Parameter(TrackParticle::dz);
// //   std::cout<<"deb 215"<<std::endl; 
//   double phi = MuTrack.Parameter(TrackParticle::phi);
// //   std::cout<<"deb 216"<<std::endl; 
//   std::cout<<"muon parameters  (phi, dz , dxy )"<< phi<< "  "<<dz<<"  "<<"  "<<dxy<<std::endl;  
//   double ddxy  = sqrt(MuTrack.Covariance(TrackParticle::dxy,TrackParticle::dxy));
// //   std::cout<<"deb 217"<<std::endl; 
//    double ddz   = sqrt(MuTrack.Covariance(TrackParticle::dz,TrackParticle::dz));
// //   std::cout<<"deb 218"<<std::endl; 
//    double ddphi = sqrt(MuTrack.Covariance(TrackParticle::phi,TrackParticle::phi));
// //   std::cout<<"deb 219"<<std::endl; 
//   double theta = asin(dxy/(sqrt(dxy*dxy + dz*dz)));
// //   std::cout<<"deb 210"<<std::endl; 
// //   std::cout<<"muon parameters theta "<< theta<<std::endl;  
// //   std::cout<<"muon parameters ddphi   ddz  ddxy "<< ddphi<< "  "<<ddz<<"  "<<"  "<<ddxy<<std::endl; 
//   double dtheta = sqrt(dz*dz*dz*dz*ddxy*ddxy/(dxy*dxy + dz*dz)/(dxy*dxy + dz*dz)/(dxy*dxy + dz*dz)   + dz*dz*dxy*dxy*ddz*ddz/(dxy*dxy + dz*dz)/(dxy*dxy + dz*dz)/(dxy*dxy + dz*dz))/sin(theta);
// //   std::cout<<"deb 2111"<<std::endl; 
//   double dx =  sqrt(cos(theta)*cos(theta)*cos(phi)*cos(phi)*dtheta*dtheta   + sin(theta)*sin(theta)*sin(phi)*sin(phi)*ddphi*ddphi );
// //   std::cout<<"deb 2112"<<std::endl; 
//   double dy =  sqrt(cos(theta)*cos(theta)*sin(phi)*sin(phi)*dtheta*dtheta   + sin(theta)*sin(theta)*cos(phi)*cos(phi)*ddphi*ddphi );
// //   std::cout<<"deb 2113"<<std::endl; 
//    double deltaz =  sin(theta)*dtheta;  
// //   std::cout<<"deb 2114"<<std::endl; 



//   TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
//   TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);
// //   std::cout<<"==========================="<<std::endl; 
//   par(LorentzVectorParticle::vx,0)=MuTrack.Parameter(LorentzVectorParticle::vx);
//   par(LorentzVectorParticle::vy,0)=MuTrack.Parameter(LorentzVectorParticle::vy);
//   par(LorentzVectorParticle::vz,0)=MuTrack.Parameter(LorentzVectorParticle::vz);
//   par(LorentzVectorParticle::px,0)=1;
//   par(LorentzVectorParticle::py,0)=1;
//   par(LorentzVectorParticle::pz,0)=1;
//   par(LorentzVectorParticle::m,0) =0.105;




//   Cov(LorentzVectorParticle::vx,LorentzVectorParticle::vx) = dx*dx;
//   Cov(LorentzVectorParticle::vy,LorentzVectorParticle::vy) = dy*dy;
//   Cov(LorentzVectorParticle::vz,LorentzVectorParticle::vz) = deltaz*deltaz;

//   Cov(LorentzVectorParticle::vx,LorentzVectorParticle::vy) = dx*dy;
//   Cov(LorentzVectorParticle::vx,LorentzVectorParticle::vz) = dx*deltaz;


//   Cov(LorentzVectorParticle::vy,LorentzVectorParticle::vx) = dy*dx;
//   Cov(LorentzVectorParticle::vy,LorentzVectorParticle::vz) = dy*deltaz;


//   Cov(LorentzVectorParticle::vz,LorentzVectorParticle::vx) = deltaz*dx;
//   Cov(LorentzVectorParticle::vz,LorentzVectorParticle::vy) = deltaz*dy;
//    std::cout<<"deb 22--->"<<std::endl; 
// //   //----
// //    std::cout<<"ConvertTrackParticleToLorentzVectorParticle   Cov Muon  "<<std::endl;
// //   for(int str =0; str < Cov.GetNrows(); str++){ 
// //     for(int kol =0; kol < Cov.GetNcols(); kol++){
// //        std::cout<<"  "<< Cov(str,kol)<<"  ";

// //     }    
// //      std::cout<<std::endl; 
// //   }

//   return LorentzVectorParticle(par,Cov,PDGInfo::mu_minus,0,MuTrack.BField());

// }


// LorentzVectorParticle DiTauConstrainedFitter::TauMuEstimatorHelixApproximation(TrackParticle MuTrack, TLorentzVector MuonLV, TVector3 PV,TVector3 SV,TVector3 TauDir,TVector3 TauDirError, LorentzVectorParticle Mu, double TauA1Energy, TMatrixTSym<double> VertexCov){
//   std::vector<double> outxyz;
  

//      double phi = TauDir.Phi();
//      double phidelta = (phi*phi +1  )*sqrt(TauDirError.Y()*TauDirError.Y()/TauDir.X()/TauDir.X()    +  TauDir.Y()*TauDir.Y()* TauDirError.X()*TauDirError.X()/TauDir.X()/TauDir.X()/TauDir.X()/TauDir.X());
//   //   double ox,oy,oz;
  
//   //   ox = TauDirError.X();
//   //   oy = TauDirError.Y();
//   //   oz = TauDirError.Z();

  
//   double dxy   =MuTrack.Parameter(TrackParticle::dxy);
//   double kappa =MuTrack.Parameter(TrackParticle::kappa);
//   double phi0  =MuTrack.Parameter(TrackParticle::phi);
//   double lam   =MuTrack.Parameter(TrackParticle::lambda);
//   double dz    =MuTrack.Parameter(TrackParticle::dz);
//   double c     =MuTrack.Charge();

//   double xpoca2 = dxy*sin(phi0);
//   double ypoca2 = -dxy*cos(phi0);


//   //   double alpha = c/5.24784503128485298e+01;

//   double rho = kappa/MuTrack.BField();
//   double a1 = -c*4*2.9998/1000;



//   double px0 = cos(phi0)/fabs(rho);
//   double py0 = sin(phi0)/fabs(rho);
//   double pz0 = tan(lam)/fabs(rho);

  
//   double xc = dxy*sin(phi0) - py0/a1, yc = dxy*cos(phi0) + px0/a1;
//   double PhiToHelixCenter = atan((yc + dxy*cos(phi0))/ ( xc - dxy*sin(phi0)));
   

//   TVector3  Po = IntersectionPoint(PV, SV, xc - PV.X(),  yc - PV.Y(),MuonLV.Pt()/a1 );
//   if(Po.X()!=0 and Po.Y()!=0){
//     double PhiAtPoint = atan((Po.Y()-yc)/(Po.X()-xc));
//     std::cout<<"PhiAtPoint "<<PhiAtPoint<<std::endl;
//     double b= MuTrack.BField();
//     double rad = 1/rho/a1;
//     double sTPo = (PhiAtPoint - PhiToHelixCenter)*rad;
//     double deltaphi = PhiAtPoint - PhiToHelixCenter;
//     double rg   =  MuonLV.Pt()/a1;

// //     double xPoint =  (dxy -rg)*sin(phi0) +  rg*sin(deltaphi + phi0);
// //     double yPoint =  (dxy +rg)*cos(phi0) -  rg*cos(deltaphi + phi0);
// //     double zPoint =  dz  +  tan(lam)*rg*deltaphi;
//     double phi1 = atan(Po.X()*tan(TauDir.Phi())  - yc)/(Po.Y() - xc);
//     //    double phi1 = atan2(Po.X(),Po.Y());

//     double xPoint =  (dxy - rg)*sin(phi0)  + rg* sin(phi1);
//     double yPoint =  (dxy+rg)*cos(phi0) - rg*cos(phi1);
//     double zPoint =  dz  +  tan(lam)*rg*deltaphi;

//     double phi2check = asin(sin(phi0) - (Po.X() - xpoca2)/c/rg);

// //     double xPoint2 =  xpoca2 + c* rg*(-sin(phi2check) + sin(phi0));
// //     double yPoint2 =  ypoca2 + c* rg*(cos(phi2check) - cos(phi0));
// //     double zPoint2 = dz - c*(phi2check -phi0 )*rg*tan(lam);


//     double xPoint2 =  xpoca2 + c* rg*(-sin(phi1) + sin(phi0));
//     double yPoint2 =  ypoca2 + c* rg*(cos(phi1) - cos(phi0));
//     double zPoint2 = dz - c*(phi1 -phi0 )*rg*tan(lam);
  

// //   std::cout<<" alternative x,y,z  at the curvature IntersectionPoint "<<xPoint2<< "  "<< yPoint2<<"   " << zPoint2<<std::endl;

// //     std::cout<<" x,y,z  at the curvature  "<<xPoint<< "  "<< yPoint<<"   " << zPoint<<std::endl;
// //     std::cout<<" alternative x,y,z  at the curvature  "<<xPoint2<< "  "<< yPoint2<<"   " << zPoint2<<std::endl;
// //     std::cout<<" compare two phis  "<<PhiAtPoint<< "  "<<phi1 <<std::endl;

//     double phifin;
//     double pi = 3.14159265359;
//     TVector3 TauMuPoint(xPoint,yPoint,zPoint);
    
//     double x2Corr = xPoint - PV.X();
//     double y2Corr = yPoint - PV.Y();
//     double z2Corr = zPoint - PV.Z();
// //     std::cout<<"PV   "<<PV.X()<<"  "<<PV.Y()<<"  "<<PV.Z()<<"   uncorrected  phi " << "  "<< atan2(yPoint,xPoint)<<"  correctred   "<<atan2(y2Corr,x2Corr)  <<"   TauMuPoint  " <<TauMuPoint.Phi() <<std::endl;
// //     std::cout<<"ReturnDecayPoint--->  DP3   "<<x2Corr<<"  "<<y2Corr<<"  "<<z2Corr<<"   uncorrected  phi " << "  "<< atan2(yPoint,xPoint)<<"  correctred   "<<atan2(y2Corr,x2Corr)  <<"   TauMuPoint  " <<TauMuPoint.Phi() <<std::endl;
// //     std::cout<<" FlkightLength      "<< sqrt(xPoint*xPoint +yPoint*yPoint+ zPoint*zPoint) <<std::endl;

//     std::cout<<" TauDir.Phi()  "<<TauDir.Phi()<<std::endl;
//     TMatrixT<double>  CovDerivative;
//     CovDerivative.ResizeTo(3,6);

//     double kap = kappa*5.24784503128485298e+01;
//     double r = 1/kap;
//     CovDerivative(0,TrackParticle::kappa) = -sin(phi0) + sin(deltaphi + phi0);
//     CovDerivative(0,TrackParticle::lambda) =0;
//     CovDerivative(0,TrackParticle::phi) = (dxy -rg)*cos(phi0) +  rg*cos(deltaphi + phi0);
//     CovDerivative(0,TrackParticle::dxy) = sin(PhiToHelixCenter);
//     CovDerivative(0,TrackParticle::dz) = 0;
//     CovDerivative(0,5) = 1;//r*cos(deltaphi + phi0);
    
//     CovDerivative(1,TrackParticle::kappa) = cos(phi0) - cos(deltaphi + phi0);
//     CovDerivative(1,TrackParticle::lambda) = 0;
//     CovDerivative(1,TrackParticle::phi) = -(r+dxy)*sin(phi0)  + r*sin(deltaphi + phi0);
//     CovDerivative(1,TrackParticle::dxy) = cos(phi0);
//     CovDerivative(1,TrackParticle::dz) = 0;
//     CovDerivative(1,5) = 1;//r*sin(deltaphi + phi0);

    
//     CovDerivative(2,TrackParticle::kappa) =lam*deltaphi;
//     CovDerivative(2,TrackParticle::lambda) = r*deltaphi;
//     CovDerivative(2,TrackParticle::phi) = 0;
//     CovDerivative(2,TrackParticle::dxy) = 0;
//     CovDerivative(2,TrackParticle::dz) =1;
//     CovDerivative(2,5) =1;//r*lam;

//     TMatrixT<double>  CovHelix;
//     CovHelix.ResizeTo(TrackParticle::NHelixPar+1,TrackParticle::NHelixPar+1);
//     for(int i=0; i<TrackParticle::NHelixPar;i++){
//       for(int j=0; j<TrackParticle::NHelixPar;j++){
// 	CovHelix(i,j) = MuTrack.Covariance(i,j);
//       }
//     }
//     CovHelix(0,5) = 0;
//     CovHelix(1,5) = 0;
//     CovHelix(2,5) = 0;
//     CovHelix(3,5) = 0;
//     CovHelix(4,5) = 0;
//     CovHelix(0,5) = 0;
    
//     CovHelix(5,0) = 0;
//     CovHelix(5,1) = 0;
//     CovHelix(5,2) = 0;
//     CovHelix(5,3) = 0;
//     CovHelix(5,4) = 0;
    

//     CovHelix(5,5) = 1;//phidelta*phidelta;


//        TMatrixT<double> CovDerivativeT=CovDerivative; CovDerivativeT.T();
//        TMatrixT<double> CovVertexXYZ  = CovDerivative*CovHelix*CovDerivativeT;
    
//        TMatrixT<double> CovVertexXYZ22  = CovHelix*CovDerivativeT;
//        TMatrixT<double> CovVertexXYZ33  = CovDerivative*CovVertexXYZ22;
    


//     //    ////----
//     std::cout<<"CovDerivative "<<std::endl;
//        for(int str =0; str < CovDerivative.GetNrows(); str++){
//        for(int kol =0; kol < CovDerivative.GetNcols(); kol++){
//             std::cout<<"  "<< CovDerivative(str,kol)<<"  ";
	    
// 	}    
// 	std::cout<<std::endl;
//      }

//         ////----
//         std::cout<<"CovDerivativeT "<<std::endl;
//        for(int str =0; str < CovDerivativeT.GetNrows(); str++){
//         for(int kol =0; kol < CovDerivativeT.GetNcols(); kol++){
//             std::cout<<"  "<< CovDerivativeT(str,kol)<<"  ";
    
//          }    
//          std::cout<<std::endl;
//        }


//     //    ////----
//     //    std::cout<<"InitialCovariance "<<std::endl;
//     //   for(int str =0; str < CovDerivativeT.GetNrows(); str++){
//     //     for(int kol =0; kol < CovDerivativeT.GetNcols(); kol++){
//     //        std::cout<<"  "<< CovDerivativeT(str,kol)<<"  ";
    
//     //     }    
//     //     std::cout<<std::endl;
//     //   }



// //         ////----
// //         std::cout<<"DiTauConstraint V_helix "<<std::endl;
// //        for(int str =0; str < CovHelix.GetNrows(); str++){
// //          for(int kol =0; kol < CovHelix.GetNcols(); kol++){
// //             std::cout<<"  "<< CovHelix(str,kol)<<"  ";
// //          }    
// //          std::cout<<std::endl;
// //        }
    
// //     //   ////----
// //     //   std::cout<<"DiTauConstraint initial par covariance "<<std::endl;
// //     //   for(int i=0; i<TrackParticle::NHelixPar;i++){
// //     //     for(int j=0; j<TrackParticle::NHelixPar;j++){
// //     //       std::cout<<"  "<< MuTrack.Covariance(i,j)<<"  ";
    
// //     //     }
// //     //    std::cout<<std::endl;
// //     //   }


// //        std::cout<<"DiTauConstraint CovVertexXYZ "<<std::endl;
// //        for(int str =0; str < CovVertexXYZ.GetNrows(); str++){
// //          for(int kol =0; kol < CovVertexXYZ.GetNcols(); kol++){
// //             std::cout<<"  "<< CovVertexXYZ(str,kol)<<"  ";
    
// //          }    
// //          std::cout<<std::endl;
// //        }
    
// //        std::cout<<"DiTauConstraint CovVertexXYZ22 "<<std::endl;
// //        for(int str =0; str < CovVertexXYZ22.GetNrows(); str++){
// //          for(int kol =0; kol < CovVertexXYZ22.GetNcols(); kol++){
// //             std::cout<<"  "<< CovVertexXYZ22(str,kol)<<"  ";
// //             }    
// //          std::cout<<std::endl;
// //        }
    
// //     //   std::cout<<"DiTauConstraint CovDerivative "<<std::endl;
// //     //   for(int str =0; str < CovDerivative.GetNrows(); str++){
// //     //     for(int kol =0; kol < CovDerivative.GetNcols(); kol++){
// //     //        std::cout<<"  "<< CovDerivative(str,kol)<<"  ";
    
// //     //     }    
// //     //     std::cout<<std::endl;
// //     //   }


// //        std::cout<<"DiTauConstraint CovVertexXYZ33 "<<std::endl;
// //        for(int str =0; str < CovVertexXYZ33.GetNrows(); str++){
// //          for(int kol =0; kol < CovVertexXYZ33.GetNcols(); kol++){
// //             std::cout<<"  "<< CovVertexXYZ33(str,kol)<<"  ";
    
// //          }    
// //          std::cout<<std::endl;
// //        }
//     TVector3 outVector;
//     TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
//     TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);
//     //   std::cout<<"==========================="<<std::endl; 
//     double rApml = sqrt(pow(xPoint,2) + pow(yPoint,2) + pow(zPoint,2));
//     double cosPhi =( TauDir.X()*xPoint + TauDir.Z()*yPoint + TauDir.Z()*zPoint)/rApml/TauDir.Mag();
    
//     TRandom *random = new TRandom();  //  Generate Random value around z mass
//     double RZMass = random->Gaus(91.18,2.5);
//     double TauMuEnergyEstimate = 0.5*RZMass*RZMass/TauA1Energy/(1 - cosPhi);
//     TLorentzVector TestTauMu(TauMuEnergyEstimate*xPoint/rApml,TauMuEnergyEstimate*yPoint/rApml,TauMuEnergyEstimate*zPoint/rApml,TauMuEnergyEstimate);

//     std::cout<<"-=-=-=-=-=-=-=-=-=-   Test pHi "<<TestTauMu.Phi() <<std::endl; 
//     par(LorentzVectorParticle::vx,0)=xPoint;
//     par(LorentzVectorParticle::vy,0)=yPoint;
//     par(LorentzVectorParticle::vz,0)=zPoint;
//     par(LorentzVectorParticle::px,0)=TauMuEnergyEstimate*xPoint/rApml;
//     par(LorentzVectorParticle::py,0)=TauMuEnergyEstimate*yPoint/rApml;
//     par(LorentzVectorParticle::pz,0)=TauMuEnergyEstimate*zPoint/rApml;
//     par(LorentzVectorParticle::m,0) =1.777;

//     for(int str =0; str < 3 ; str++){ 
//       for(int kol =0; kol < 3; kol++){
// 	Cov(str,kol) = CovVertexXYZ(str,kol);
	
//       }    
//     }
    
//     for(int str =3; str < 6; str++){ 
//       for(int kol =3; kol < 6; kol++){
// 	Cov(str,kol) = CovVertexXYZ(str,kol);
	
//       }    
//     }
//     std::cout<<"COV "<<std::endl;
//     for(int str =0; str < Cov.GetNrows(); str++){
//       for(int kol =0; kol < Cov.GetNcols(); kol++){
// 	std::cout<<"  "<< Cov(str,kol)<<"  ";
	
//       }    
//       std::cout<<std::endl;
//     }
    
    

//     return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);
    
//   }
// }
 
// LorentzVectorParticle DiTauConstrainedFitter::TauMuEstimatorLinearApproximation(TrackParticle MuTrack, TLorentzVector MuonLV, TVector3 PV,TVector3 SV,TVector3 TauDir,TVector3 TauDirError, LorentzVectorParticle Mu, double TauA1Energy, TMatrixTSym<double> VertexCov,TVector3 MuonPoca,LorentzVectorParticle TauA1){
//    std::vector<double> outxyz;

//       double dxy   =MuTrack.Parameter(TrackParticle::dxy);
//       double kappa =MuTrack.Parameter(TrackParticle::kappa);
//       double phi0  =MuTrack.Parameter(TrackParticle::phi);
//       double lam   =MuTrack.Parameter(TrackParticle::lambda);
//       double dz    =MuTrack.Parameter(TrackParticle::dz);
//       double c     =MuTrack.Charge();


//       double phi = TauDir.Phi();
//       double phidelta = (phi*phi +1  )*sqrt(TauDirError.Y()*TauDirError.Y()/TauDir.X()/TauDir.X()    +  TauDir.Y()*TauDir.Y()* TauDirError.X()*TauDirError.X()/TauDir.X()/TauDir.X()/TauDir.X()/TauDir.X());
//    //   double ox,oy,oz;

//    //   ox = TauDirError.X();
//    //   oy = TauDirError.Y();
//    //   oz = TauDirError.Z();


//       TVector3 ReferencePoint  = MuonPoca;

//       double aMuonTrack  = atan2(ReferencePoint.Y(), ReferencePoint.X());
//       double bMuonTrack = ReferencePoint.Y() - aMuonTrack*ReferencePoint.X();


//    //   double alpha = c/5.24784503128485298e+01;

//       double rho = kappa/MuTrack.BField();
//       double a1 = -c*4*2.9998/1000;



//    double px0 = cos(phi0)/fabs(rho);
//    double py0 = sin(phi0)/fabs(rho);
//    double pz0 = tan(lam)/fabs(rho);


//    double xc = dxy*sin(phi0) - py0/a1, yc = dxy*cos(phi0) + px0/a1;
//    double PhiToHelixCenter = atan((yc - dxy*cos(phi0))/ ( xc - dxy*sin(phi0)));
 
//    std::cout<<"Muon paramters "<< c*dxy*sin(phi0) << "  "<< -c*dxy*cos(phi0)<<std::endl;
//    TVector3  Po(1,1,1);// = IntersectionPointLinearApproximation( PV, SV,  MuonPoca,  MuonLV, MuTrack,TauDirError,TauA1);
//    if(Po.X()!=0 and Po.Y()!=0){
//      double PhiAtPoint = atan((Po.Y()-yc)/(Po.X()-xc));
//      std::cout<<"PhiAtPoint "<<PhiAtPoint<<std::endl;


//      double b= MuTrack.BField();
//      double rad = 1/rho/a1;
//      double sTPo = (PhiAtPoint - PhiToHelixCenter)*rad;
//      double deltaphi = PhiAtPoint - PhiToHelixCenter;
//      double rg   =  MuonLV.Pt()/a1;

//      double xPoint =  (dxy -rg)*sin(phi0) +  rg*sin(deltaphi + phi0);
//      double yPoint =  (dxy +rg)*cos(phi0) -  rg*cos(deltaphi + phi0);
//      double zPoint =  dz  +  tan(lam)*rg*deltaphi;



//      std::cout<<" x,y,z  at the curvature  "<<xPoint<< "  "<< yPoint<<"   " << zPoint<<std::endl;


//      double phifin;
//      double pi = 3.14159265359;
//      TVector3 TauMuPoint(xPoint,yPoint,zPoint);
  
//      double x2Corr = xPoint - PV.X();
//      double y2Corr = yPoint - PV.Y();
//      double z2Corr = zPoint - PV.Z();
//  //     std::cout<<"PV   "<<PV.X()<<"  "<<PV.Y()<<"  "<<PV.Z()<<"   uncorrected  phi " << "  "<< atan2(yPoint,xPoint)<<"  correctred   "<<atan2(y2Corr,x2Corr)  <<"   TauMuPoint  " <<TauMuPoint.Phi() <<std::endl;
//  //     std::cout<<"ReturnDecayPoint--->  DP3   "<<x2Corr<<"  "<<y2Corr<<"  "<<z2Corr<<"   uncorrected  phi " << "  "<< atan2(yPoint,xPoint)<<"  correctred   "<<atan2(y2Corr,x2Corr)  <<"   TauMuPoint  " <<TauMuPoint.Phi() <<std::endl;
//  //     std::cout<<" FlkightLength      "<< sqrt(xPoint*xPoint +yPoint*yPoint+ zPoint*zPoint) <<std::endl;

//      std::cout<<" TauDir.Phi()  "<<TauDir.Phi()<<std::endl;
//      TMatrixT<double>  CovDerivative;
//      CovDerivative.ResizeTo(3,6);

//      double kap = kappa*5.24784503128485298e+01;
//      double r = 1/kap;
//      CovDerivative(0,TrackParticle::kappa) = -sin(phi0) + sin(deltaphi + phi0);
//      CovDerivative(0,TrackParticle::lambda) =0;
//      CovDerivative(0,TrackParticle::phi) = (dxy -rg)*cos(phi0) +  rg*cos(deltaphi + phi0);
//      CovDerivative(0,TrackParticle::dxy) = sin(PhiToHelixCenter);
//      CovDerivative(0,TrackParticle::dz) = 0;
//      CovDerivative(0,5) = 1;//r*cos(deltaphi + phi0);
  
//      CovDerivative(1,TrackParticle::kappa) = cos(phi0) - cos(deltaphi + phi0);
//      CovDerivative(1,TrackParticle::lambda) = 0;
//      CovDerivative(1,TrackParticle::phi) = -(r+dxy)*sin(phi0)  + r*sin(deltaphi + phi0);
//      CovDerivative(1,TrackParticle::dxy) = cos(phi0);
//      CovDerivative(1,TrackParticle::dz) = 0;
//      CovDerivative(1,5) = 1;//r*sin(deltaphi + phi0);

  
//      CovDerivative(2,TrackParticle::kappa) =lam*deltaphi;
//      CovDerivative(2,TrackParticle::lambda) = r*deltaphi;
//      CovDerivative(2,TrackParticle::phi) = 0;
//      CovDerivative(2,TrackParticle::dxy) = 0;
//      CovDerivative(2,TrackParticle::dz) =1;
//      CovDerivative(2,5) =1;//r*lam;

//      TMatrixT<double>  CovHelix;
//      CovHelix.ResizeTo(TrackParticle::NHelixPar+1,TrackParticle::NHelixPar+1);
//      for(int i=0; i<TrackParticle::NHelixPar;i++){
//        for(int j=0; j<TrackParticle::NHelixPar;j++){
//  	CovHelix(i,j) = MuTrack.Covariance(i,j);
//        }
//      }
//      CovHelix(0,5) = 0;
//      CovHelix(1,5) = 0;
//      CovHelix(2,5) = 0;
//      CovHelix(3,5) = 0;
//      CovHelix(4,5) = 0;
//      CovHelix(0,5) = 0;
  
//      CovHelix(5,0) = 0;
//      CovHelix(5,1) = 0;
//      CovHelix(5,2) = 0;
//      CovHelix(5,3) = 0;
//      CovHelix(5,4) = 0;
  

//      CovHelix(5,5) = 1;//phidelta*phidelta;


//         TMatrixT<double> CovDerivativeT=CovDerivative; CovDerivativeT.T();
//         TMatrixT<double> CovVertexXYZ  = CovDerivative*CovHelix*CovDerivativeT;
  
//         TMatrixT<double> CovVertexXYZ22  = CovHelix*CovDerivativeT;
//         TMatrixT<double> CovVertexXYZ33  = CovDerivative*CovVertexXYZ22;
  


// //      std::cout<<"TauA1.Covariance "<<std::endl;
// //         for(int str =0; str < CovDerivative.GetNrows(); str++){
// // 	  for(int kol =0; kol < TauA1.Covariance().GetNcols(); kol++){
// //              std::cout<<"  "<< TauA1.Covariance(str,kol)<<"  ";
	    
// //  	}    
// //  	std::cout<<std::endl;
// //       }


//      //    ////----
// //      std::cout<<"CovDerivative "<<std::endl;
// //         for(int str =0; str < CovDerivative.GetNrows(); str++){
// //         for(int kol =0; kol < CovDerivative.GetNcols(); kol++){
// //              std::cout<<"  "<< CovDerivative(str,kol)<<"  ";
	    
// //  	}    
// //  	std::cout<<std::endl;
// //       }


// //          ////----
// //          std::cout<<"CovDerivativeT "<<std::endl;
// //         for(int str =0; str < CovDerivativeT.GetNrows(); str++){
// //          for(int kol =0; kol < CovDerivativeT.GetNcols(); kol++){
// //              std::cout<<"  "<< CovDerivativeT(str,kol)<<"  ";
  
// //           }    
// //           std::cout<<std::endl;
// //         }


//      //    ////----
//      //    std::cout<<"InitialCovariance "<<std::endl;
//      //   for(int str =0; str < CovDerivativeT.GetNrows(); str++){
//      //     for(int kol =0; kol < CovDerivativeT.GetNcols(); kol++){
//      //        std::cout<<"  "<< CovDerivativeT(str,kol)<<"  ";
//      //     }    
//      //     std::cout<<std::endl;
//      //   }



// //          ////----
// //          std::cout<<"DiTauConstraint V_helix "<<std::endl;
// //         for(int str =0; str < CovHelix.GetNrows(); str++){
// //           for(int kol =0; kol < CovHelix.GetNcols(); kol++){
// //              std::cout<<"  "<< CovHelix(str,kol)<<"  ";
// //           }    
// //           std::cout<<std::endl;
// //         }
  
//      //   ////----
//      //   std::cout<<"DiTauConstraint initial par covariance "<<std::endl;
//      //   for(int i=0; i<TrackParticle::NHelixPar;i++){
//      //     for(int j=0; j<TrackParticle::NHelixPar;j++){
//      //       std::cout<<"  "<< MuTrack.Covariance(i,j)<<"  ";
  
//      //     }
//      //    std::cout<<std::endl;
//      //   }


// //         std::cout<<"DiTauConstraint CovVertexXYZ "<<std::endl;
// //         for(int str =0; str < CovVertexXYZ.GetNrows(); str++){
// //           for(int kol =0; kol < CovVertexXYZ.GetNcols(); kol++){
// //              std::cout<<"  "<< CovVertexXYZ(str,kol)<<"  ";
  
// //           }    
// //           std::cout<<std::endl;
// //         }
  
// //         std::cout<<"DiTauConstraint CovVertexXYZ22 "<<std::endl;
// //         for(int str =0; str < CovVertexXYZ22.GetNrows(); str++){
// //           for(int kol =0; kol < CovVertexXYZ22.GetNcols(); kol++){
// //              std::cout<<"  "<< CovVertexXYZ22(str,kol)<<"  ";
// //              }    
// //           std::cout<<std::endl;
// //         }
  
//      //   std::cout<<"DiTauConstraint CovDerivative "<<std::endl;
//      //   for(int str =0; str < CovDerivative.GetNrows(); str++){
//      //     for(int kol =0; kol < CovDerivative.GetNcols(); kol++){
//      //        std::cout<<"  "<< CovDerivative(str,kol)<<"  ";
  
//      //     }    
//      //     std::cout<<std::endl;
//      //   }


// //         std::cout<<"DiTauConstraint CovVertexXYZ33 "<<std::endl;
// //         for(int str =0; str < CovVertexXYZ33.GetNrows(); str++){
// //           for(int kol =0; kol < CovVertexXYZ33.GetNcols(); kol++){
// //              std::cout<<"  "<< CovVertexXYZ33(str,kol)<<"  ";
  
// //           }    
// //           std::cout<<std::endl;
// //         }
  


//      TVector3 outVector;



//    TMatrixT<double>    par(LorentzVectorParticle::NLorentzandVertexPar,10);
//    TMatrixTSym<double> Cov(LorentzVectorParticle::NLorentzandVertexPar);
//  //   std::cout<<"==========================="<<std::endl; 
//    double rApml = sqrt(pow(xPoint,2) + pow(yPoint,2) + pow(zPoint,2));
//    double cosPhi =( TauDir.X()*xPoint + TauDir.Z()*yPoint + TauDir.Z()*zPoint)/rApml/TauDir.Mag();

//    TRandom *random = new TRandom();  //  Generate Random value around z mass
//    double RZMass = random->Gaus(91.18,2.5);
//    double TauMuEnergyEstimate = 0.5*RZMass*RZMass/TauA1Energy/(1 - cosPhi);

//    par(LorentzVectorParticle::vx,0)=xPoint;
//    par(LorentzVectorParticle::vy,0)=yPoint;
//    par(LorentzVectorParticle::vz,0)=zPoint;
//    par(LorentzVectorParticle::px,0)=-TauA1.LV().Px();//TauMuEnergyEstimate*xPoint/rApml;
//    par(LorentzVectorParticle::py,0)=-TauA1.LV().Py();//TauMuEnergyEstimate*yPoint/rApml;
//    par(LorentzVectorParticle::pz,0)=-TauA1.LV().Pz();//TauMuEnergyEstimate*zPoint/rApml;
//    par(LorentzVectorParticle::m,0) =1.777;

//    for(int str =0; str < 3 ; str++){ 
//      for(int kol =0; kol < 3; kol++){
//        Cov(str,kol) = CovVertexXYZ(str,kol);

//      }    
//    }

//    for(int str =3; str < 6; str++){ 
//      for(int kol =3; kol < 6; kol++){
//        Cov(str,kol) = CovVertexXYZ(str,kol);
    
//      }    
//    }
//         std::cout<<"COV "<<std::endl;
//         for(int str =0; str < Cov.GetNrows(); str++){
//           for(int kol =0; kol < Cov.GetNcols(); kol++){
//              std::cout<<"  "<< Cov(str,kol)<<"  ";
  
//           }    
//           std::cout<<std::endl;
//         }
  


//   return LorentzVectorParticle(par,Cov,PDGInfo::tau_minus,0,0);

//    }
//  }
 


// std::vector<double>  
// DiTauConstrainedFitter::ReturnDecayPoint(TrackParticle MuTrack, TLorentzVector MuonLV, TVector3 PV,TVector3 SV,TVector3 TauDir,TVector3 TauDirError){
//   std::vector<double> outxyz;
  

//      double phi = TauDir.Phi();
//      double phidelta = (phi*phi +1  )*sqrt(TauDirError.Y()*TauDirError.Y()/TauDir.X()/TauDir.X()    +  TauDir.Y()*TauDir.Y()* TauDirError.X()*TauDirError.X()/TauDir.X()/TauDir.X()/TauDir.X()/TauDir.X());
//   //   double ox,oy,oz;
  
//   //   ox = TauDirError.X();
//   //   oy = TauDirError.Y();
//   //   oz = TauDirError.Z();

  
//   double dxy   =MuTrack.Parameter(TrackParticle::dxy);
//   double kappa =MuTrack.Parameter(TrackParticle::kappa);
//   double phi0  =MuTrack.Parameter(TrackParticle::phi);
//   double lam   =MuTrack.Parameter(TrackParticle::lambda);
//   double dz    =MuTrack.Parameter(TrackParticle::dz);
//   double c     =MuTrack.Charge();

//   //   double alpha = c/5.24784503128485298e+01;
//   //  std::cout<<"  "
//   double rho = kappa/MuTrack.BField();
//   double a1 = -c*4*2.9998/1000;

//   double px0 = cos(phi0)/fabs(rho);
//   double py0 = sin(phi0)/fabs(rho);
//   double pz0 = tan(lam)/fabs(rho);

  
//   double xc = dxy*sin(phi0) - py0/a1, yc = dxy*cos(phi0) + px0/a1;
//   double PhiToHelixCenter = atan((yc - dxy*cos(phi0))/ ( xc - dxy*sin(phi0)));
   

//   TVector3  Po = IntersectionPoint(PV, SV, xc - PV.X(),  yc - PV.Y(),MuonLV.Pt()/a1 );
//   if(Po.X()!=0 and Po.Y()!=0){
//     double PhiAtPoint = atan((Po.Y()-yc)/(Po.X()-xc));
//     std::cout<<"PhiAtPoint "<<PhiAtPoint<<std::endl;
//     double b= MuTrack.BField();
//     double rad = 1/rho/a1;
//     double sTPo = (PhiAtPoint - PhiToHelixCenter)*rad;
//     double deltaphi = PhiAtPoint - PhiToHelixCenter;
//     double rg   =  MuonLV.Pt()/a1;

//     double xPoint =  (dxy -rg)*sin(phi0) +  rg*sin(deltaphi + phi0);
//     double yPoint =  (dxy +rg)*cos(phi0) -  rg*cos(deltaphi + phi0);
//     double zPoint =  dz  +  tan(lam)*rg*deltaphi;

//     // std::cout<<" x,y,z  at the curvature  "<<xPoint<< "  "<< yPoint<<"   " << zPoint<<std::endl;

//     double phifin;
//     double pi = 3.14159265359;
//     TVector3 TauMuPoint(xPoint,yPoint,zPoint);
    
//     double x2Corr = xPoint - PV.X();
//     double y2Corr = yPoint - PV.Y();
//     double z2Corr = zPoint - PV.Z();

//     //std::cout<<"PV   "<<PV.X()<<"  "<<PV.Y()<<"  "<<PV.Z()<<"   uncorrected  phi " << "  "<< atan2(yPoint,xPoint)<<"  correctred   "<<atan2(y2Corr,x2Corr)  <<"   TauMuPoint  " <<TauMuPoint.Phi() <<std::endl;
//     //std::cout<<"ReturnDecayPoint--->  DP3   "<<x2Corr<<"  "<<y2Corr<<"  "<<z2Corr<<"   uncorrected  phi " << "  "<< atan2(yPoint,xPoint)<<"  correctred   "<<atan2(y2Corr,x2Corr)  <<"   TauMuPoint  " <<TauMuPoint.Phi() <<std::endl;
//     //std::cout<<" FlkightLength      "<< sqrt(xPoint*xPoint +yPoint*yPoint+ zPoint*zPoint) <<std::endl;

//     std::cout<<" TauDir.Phi()  "<<TauDir.Phi()<<std::endl;
//     TMatrixT<double>  CovDerivative;
//     CovDerivative.ResizeTo(3,6);

//     double kap = kappa*5.24784503128485298e+01;
//     double r = 1/kap;
//     CovDerivative(0,TrackParticle::kappa) = -sin(phi0) + sin(deltaphi + phi0);
//     CovDerivative(0,TrackParticle::lambda) =0;
//     CovDerivative(0,TrackParticle::phi) = (dxy -rg)*cos(phi0) +  rg*cos(deltaphi + phi0);
//     CovDerivative(0,TrackParticle::dxy) = sin(PhiToHelixCenter);
//     CovDerivative(0,TrackParticle::dz) = 0;
//     CovDerivative(0,5) = 1;//r*cos(deltaphi + phi0);
    
//     CovDerivative(1,TrackParticle::kappa) = cos(phi0) - cos(deltaphi + phi0);
//     CovDerivative(1,TrackParticle::lambda) = 0;
//     CovDerivative(1,TrackParticle::phi) = -(r+dxy)*sin(phi0)  + r*sin(deltaphi + phi0);
//     CovDerivative(1,TrackParticle::dxy) = cos(phi0);
//     CovDerivative(1,TrackParticle::dz) = 0;
//     CovDerivative(1,5) = 1;//r*sin(deltaphi + phi0);

    
//     CovDerivative(2,TrackParticle::kappa) =lam*deltaphi;
//     CovDerivative(2,TrackParticle::lambda) = r*deltaphi;
//     CovDerivative(2,TrackParticle::phi) = 0;
//     CovDerivative(2,TrackParticle::dxy) = 0;
//     CovDerivative(2,TrackParticle::dz) =1;
//     CovDerivative(2,5) =1;//r*lam;

//        TMatrixT<double>  CovHelix;
//        CovHelix.ResizeTo(TrackParticle::NHelixPar+1,TrackParticle::NHelixPar+1);




//      for(int i=0; i<TrackParticle::NHelixPar;i++){
//          for(int j=0; j<TrackParticle::NHelixPar;j++){
//            CovHelix(i,j) = MuTrack.Covariance(i,j);
//          }
//        }
//        CovHelix(0,5) = 0;
//        CovHelix(1,5) = 0;
//        CovHelix(2,5) = 0;
//        CovHelix(3,5) = 0;
//        CovHelix(4,5) = 0;
//        CovHelix(0,5) = 0;

//        CovHelix(5,0) = 0;
//        CovHelix(5,1) = 0;
//        CovHelix(5,2) = 0;
//        CovHelix(5,3) = 0;
//        CovHelix(5,4) = 0;


//        CovHelix(5,5) = 1;//phidelta*phidelta;


//        TMatrixT<double> CovDerivativeT=CovDerivative; CovDerivativeT.T();
//        TMatrixT<double> CovVertexXYZ  = CovDerivative*CovHelix*CovDerivativeT;
    
//        TMatrixT<double> CovVertexXYZ22  = CovHelix*CovDerivativeT;
//        TMatrixT<double> CovVertexXYZ33  = CovDerivative*CovVertexXYZ22;
    


//     //    ////----
// //     std::cout<<"CovDerivative "<<std::endl;
// //        for(int str =0; str < CovDerivative.GetNrows(); str++){
// //        for(int kol =0; kol < CovDerivative.GetNcols(); kol++){
// //             std::cout<<"  "<< CovDerivative(str,kol)<<"  ";
	    
// // 	}    
// // 	std::cout<<std::endl;
// //      }

//         ////----
// //         std::cout<<"CovDerivativeT "<<std::endl;
// //        for(int str =0; str < CovDerivativeT.GetNrows(); str++){
// //         for(int kol =0; kol < CovDerivativeT.GetNcols(); kol++){
// //             std::cout<<"  "<< CovDerivativeT(str,kol)<<"  ";
    
// //          }    
// //          std::cout<<std::endl;
// //        }


//     //    ////----
//     //    std::cout<<"InitialCovariance "<<std::endl;
//     //   for(int str =0; str < CovDerivativeT.GetNrows(); str++){
//     //     for(int kol =0; kol < CovDerivativeT.GetNcols(); kol++){
//     //        std::cout<<"  "<< CovDerivativeT(str,kol)<<"  ";
    
//     //     }    
//     //     std::cout<<std::endl;
//     //   }



//         ////----
// //         std::cout<<"DiTauConstraint V_helix "<<std::endl;
// //        for(int str =0; str < CovHelix.GetNrows(); str++){
// //          for(int kol =0; kol < CovHelix.GetNcols(); kol++){
// //             std::cout<<"  "<< CovHelix(str,kol)<<"  ";
// //          }    
// //          std::cout<<std::endl;
// //        }
    
//     //   ////----
//     //   std::cout<<"DiTauConstraint initial par covariance "<<std::endl;
//     //   for(int i=0; i<TrackParticle::NHelixPar;i++){
//     //     for(int j=0; j<TrackParticle::NHelixPar;j++){
//     //       std::cout<<"  "<< MuTrack.Covariance(i,j)<<"  ";
    
//     //     }
//     //    std::cout<<std::endl;
//     //   }


// //        std::cout<<"DiTauConstraint CovVertexXYZ "<<std::endl;
// //        for(int str =0; str < CovVertexXYZ.GetNrows(); str++){
// //          for(int kol =0; kol < CovVertexXYZ.GetNcols(); kol++){
// //             std::cout<<"  "<< CovVertexXYZ(str,kol)<<"  ";
    
// //          }    
// //          std::cout<<std::endl;
// //        }
    
// //        std::cout<<"DiTauConstraint CovVertexXYZ22 "<<std::endl;
// //        for(int str =0; str < CovVertexXYZ22.GetNrows(); str++){
// //          for(int kol =0; kol < CovVertexXYZ22.GetNcols(); kol++){
// //             std::cout<<"  "<< CovVertexXYZ22(str,kol)<<"  ";
// //             }    
// //          std::cout<<std::endl;
// //        }
    
//     //   std::cout<<"DiTauConstraint CovDerivative "<<std::endl;
//     //   for(int str =0; str < CovDerivative.GetNrows(); str++){
//     //     for(int kol =0; kol < CovDerivative.GetNcols(); kol++){
//     //        std::cout<<"  "<< CovDerivative(str,kol)<<"  ";
    
//     //     }    
//     //     std::cout<<std::endl;
//     //   }
// //        std::cout<<"DiTauConstraint CovVertexXYZ33 "<<std::endl;
// //        for(int str =0; str < CovVertexXYZ33.GetNrows(); str++){
// //          for(int kol =0; kol < CovVertexXYZ33.GetNcols(); kol++){
// //             std::cout<<"  "<< CovVertexXYZ33(str,kol)<<"  ";
    
// //          }    
// //          std::cout<<std::endl;
// //        }
    


//     TVector3 outVector;
//   }
//   return outxyz;
// }



// TVector3
// DiTauConstrainedFitter::IntersectionPoint(TVector3 PV,TVector3 SV,double xc, double yc, double r){
//   TVector3 Point;
//   double xout=0,yout=0;

//   TVector3 A1TauDirection(SV.X() - PV.X(),SV.Y() - PV.Y(),SV.Z() - PV.Z()); // With respect to PV

//    std::cout<<"PhiOfTau "<< A1TauDirection.Phi()<<std::endl;
// //   std::cout<<"opposite to tau "<< tan(A1TauDirection.Phi() - 3.1415)<< "  "<<tan(A1TauDirection.Phi() + 3.1415)<<std::endl;

//   double a = tan(A1TauDirection.Phi()); 

//   double aTau = (SV.Y() - PV.Y())/(SV.X() - PV.X());
//   double bTau = (SV.X()*PV.Y() - PV.X()*SV.Y())/(SV.X() - PV.X());

// //   double A = 1 + a*a;
// //   double B = - 2*(xc + yc*a);
// //   double C = xc*xc + yc*yc - r*r;

//   double A=1+ aTau*aTau;
//   double B=2*aTau*bTau - 2*yc*aTau  - 2*xc;
//   double C=yc*yc + xc*xc - r*r + bTau*bTau -2*yc*bTau;



//   std::cout<<" A,B,C "<<A<<"  "<<B<<"  "  <<C<<std::endl; 
//   std::cout<<"D "<<B*B<<"  -  "<<4*A*C<<std::endl; 

//   double x1  = (-B + sqrt(B*B - 4*A*C))/(2*A);
//   double x2  = (-B - sqrt(B*B - 4*A*C))/(2*A);

//   double y1  = aTau*x1+bTau;
//   double y2  = aTau*x2+bTau;

//   //   std::cout<<" DiTauConstrainedFitter::IntersectionPoint  x1, y1  "<<x1<<"  "<< y1<<std::endl; 
//   //   std::cout<<" DiTauConstrainedFitter::IntersectionPoint  x2, y2  "<<x2<<"  "<< y2<<std::endl; 
 
//   if(sqrt(pow(x1,2) + pow(y1,2)) <  sqrt(pow(x2,2) + pow(y2,2))){xout = x1; yout = y1;}
//   else {xout = x2; yout = y2;}



  


//   //   std::cout<<" atan2  x2, y2  "<<atan2(yout ,xout)<<"   "<<yout <<"  =  "<<xout*a <<std::endl; 
//   //   std::cout<<" PV "<< PV.X()<<"   "<< PV.Y()<<"  SV   "<< SV.X()<<"   "<<  SV.Y() <<std::endl; 
//   //   std::cout<<" xc "<< xc<<"   "<< yc<<std::endl; 
  
//   //   std::cout<<" x2, y2  "<<xout <<"   "<<yout <<"  =  "<< SV.X() - PV.X()<<"   "<<  SV.Y()  - PV.Y()<<std::endl; 
  
//   //   double phi = atan2(yout - yc, xout - xc);
//   //   double sT = (phi - phi0)/a1/rho;
  
//   //    double xcu = x0 + (px0/a1)*sin(phi) -  (py0/a1)*(1 - cos(phi));
//   //    double ycu = y0 + (px0/a1)*(1 - cos(phi)) + (py0/a1)*sin(phi);

//   // //   double ycu = r*sin(phi) - (r + d0)*sin(phi);
//   // //   double xcu = -r*cos(phi) + (r + d0)*cos(phi);
//   //   //  std::cout<<" xn,yn at the curvature  "<<xn << "  "<< yn<<"   " <<std::endl;


//   //   std::cout<<" xcu, ycu phi  "<< xcu    <<"  "<< ycu   <<" and phi  "<< phi <<std::endl; 
//   //   std::cout<<"  "<<  ycu  <<"  =   "<< xcu*a    <<"  +  "<< b    <<"  =   "<< xcu*a + b <<std::endl; 

//    std::cout<<"comparephis helix "<<  atan2(SV.Y()-PV.Y(),SV.X()-PV.X()) <<"  == "<<  atan2(SV.Y()- yout,SV.X()- xout )<<"  "<<  "fabs "<< fabs(atan2(SV.Y()-PV.Y(),SV.X()-PV.X())-atan2(SV.Y()- yout,SV.X()-  xout))<<std::endl;
//    std::cout<<"comparephis helix "<<  yout  <<"  == "<< aTau*xout + bTau <<"  "<<  r*r <<" ==  "<< pow(yout - yc,2) + pow(xout - xc,2)<<std::endl;


//   //   std::cout<<" check both points at the curvature  "<<  sqrt(pow(x2-xc,2) + pow(y2-yc,2))   <<"  "<< sqrt(pow(x1-xc,2) + pow(y1-yc,2))   <<" and phi  "<<atan2(yout-yc,xout-xc) <<std::endl; 
  
//   //   if(B*B - 4*A*C< 0){ xout =0; yout =0;}
  
//   std::cout<<" DiTauConstrainedFitter::IntersectionPoint  xout,  yout  "<<xout<<"  "<<yout <<std::endl; 
  
//   Point.SetX(xout);
//   Point.SetY(yout);
//   Point.SetZ(0);
//   return Point;
// }





