// c++ -o read_wpluslo `root-config --glibs --cflags` CfgParser.cc LHEF.cc -lm read_test.cpp
#include "LHEF.h"
//#include "LHEF_joint.h"
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <algorithm>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "CfgParser.h"

using namespace std ;


// all histograms of a single sample
struct histos
{
  TString m_name ;  // name of the sample
  double  m_XS ;    // cross-section of the sample (in fb)
  double  m_lumi ;  // integrated luminosity (in /fb)
  map<string, TH1F *> m_histograms ;

  // integral of the events collected, 
  // needed for the histograms normalisation
  double m_norm ;
  
  histos (TString name, double XS, double lumi = 1) : 
    m_name (name), 
    m_XS (XS),
    m_lumi (lumi),
    m_norm (0.)
    {
    makeHisto ("no_shape","count","", 1,0,1) ;
    float tmp1[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
    makeHisto ("m_met","p_{t}^{miss}[GeV]","", 10,tmp1) ;
    float tmp2[]={0,20,30,40,50,60,70,80,90,100};
    makeHisto ("m_met_2","p_{t}^{miss}[GeV]","", 5,tmp2) ;
    float tmp3[]={500,1000,1500,2000,2700,3500,5000};
    float tmp6[]={0,20,40,60,80,100,120,140,160};
    float tmp9[]={5,10,20,50,80,120,160,200,250,300,370,450};
    makeHisto ("m_ptl1","p_{t}^{l1}[GeV]","", 9,tmp9) ;
    float tmp10[]={5,10,20,35,50,70,90,120,150,180,210,250};
    makeHisto ("m_ptl2","p_{t}^{l2}[GeV]","", 9,tmp10) ;
    makeHisto ("etal1","#eta_{l1}","", 20,-2.5,2.5) ;
    makeHisto ("etal2","#eta_{l2}","", 20,-2.5,2.5) ;
    }
  double increaseNorm (double step = 1.) 
    {
      m_norm += step ;
      return m_norm ;
    }

  // simplify histogram creation
  TH1F * makeHisto (const TString & varname, const TString & x_axis,const TString & y_axis, int nBins, float min, float max)
    {
      TH1F * histo = new TH1F (varname + TString ("_") + m_name, varname + TString ("_")  + m_name+TString(";")+x_axis+TString(";")+y_axis, nBins, min, max) ;
      histo->Sumw2 () ;
      if (m_histograms.find (varname.Data ()) != m_histograms.end ())
        {
          cout << "WARNING histogram " << varname << " defined twice" << endl ;
          return histo ;
        }
      m_histograms[varname.Data ()] = histo ;
      return histo ;
    } 
  TH1F * makeHisto (const TString & varname, const TString & x_axis,const TString & y_axis, int nBins, float *xbins)
    {
      TH1F * histo = new TH1F (varname + TString ("_") + m_name, varname + TString ("_")  + m_name+TString(";")+x_axis+TString(";")+y_axis, nBins, xbins) ;
      histo->Sumw2 () ;
      if (m_histograms.find (varname.Data ()) != m_histograms.end ())
        {
          cout << "WARNING histogram " << varname << " defined twice" << endl ;
          return histo ;
        }
      m_histograms[varname.Data ()] = histo ;
      return histo ;
    } 
  // fill the histograms through the map
  void fill (string varname, double value, double weight = 1.)
    { 
      bool fold_overflow=true;
      if (m_histograms.find (varname) == m_histograms.end ())
        {
          cout << "WARNING histogram " << varname << " does not exist" << endl ;
          return ;
        }
        
      int NbinsX=m_histograms[varname]->GetNbinsX();
      float BinWidth=m_histograms[varname]->GetXaxis()->GetBinWidth(NbinsX);
      float upper_edge=m_histograms[varname]->GetXaxis()->GetBinLowEdge(NbinsX)+BinWidth;
      if(value>=upper_edge && fold_overflow){
        value=m_histograms[varname]->GetXaxis()->GetBinCenter(NbinsX);
      }
      m_histograms[varname]->Fill (value, weight) ;
      return ; 
    }
  
  // normalise histograms to the integrated cross-section
  void norm (/* double inputIntegral = 0*/)
    {
//      double factor = m_lumi * m_XS / fabs (m_histograms.begin ()->second->Integral ()) ;
//      if (inputIntegral != 0) factor = m_lumi * m_XS / inputIntegral ;

      double factor = m_lumi * 1000. * m_XS / m_norm ;
      for (map<string, TH1F *>::iterator it = m_histograms.begin () ; 
           it != m_histograms.end () ; ++it)
         it->second->Scale (factor) ;
    }
  
  ~histos ()
    {
      // this being empty is not very nice, to be replaced with the unique_ptr or auto_ptr thing 
      // for the histogram pointers, but I need to find a smart way to do it
      // w/o need for implementing it for each histogram I add      
    }  
    
  void save (TFile & outfile) 
    {
      outfile.cd () ;
      for (map<string, TH1F *>::iterator it = m_histograms.begin () ; 
           it != m_histograms.end () ; ++it)
        it->second->Write () ;
    }
  
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TLorentzVector buildLP (LHEF::Reader & reader, int iPart)
{
  TLorentzVector tlv
    (
      reader.hepeup.PUP.at (iPart).at (0), //PG px
      reader.hepeup.PUP.at (iPart).at (1), //PG py
      reader.hepeup.PUP.at (iPart).at (2), //PG pz
      reader.hepeup.PUP.at (iPart).at (3)  //PG E
    ) ;
  return tlv ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


inline float zetaStar (float etaq1, float etaq2, float eta)
{
  return (eta - 0.5 * (etaq1 + etaq2)) / fabs (etaq1 - etaq2) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


// Fill the histograms for a single sample
// histograms will not get normalised, since the same sample
// could be split in several LHE files and this function
// may get called several times, one for each LHE file.
// Therefore, the normalisation will have to be called afterwards
double 
fillHistos (LHEF::Reader & reader, histos & Histos, int max = -1)
{
  int events = 0 ;
   
  //PG loop over input events
  while (reader.readEvent ()) 
    {
//      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      if (events++ % 10000 == 0) cout << "        reading event in file: " << events << endl ;
          
      vector<TLorentzVector> v_f_Ws_e ;
      vector<TLorentzVector> v_f_Ws_mu ;
      vector<TLorentzVector> v_f_quarks ;
      vector<TLorentzVector> v_f_leptons_e ;
      vector<TLorentzVector> v_f_leptons_mu ;
      vector<TLorentzVector> v_f_photons ;
      vector<TLorentzVector> v_f_neutrinos ;
      vector<int> idx_Ws_e;
      vector<int> idx_Ws_mu;
      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
          // outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
              TLorentzVector dummy = buildLP (reader, iPart) ;
              // quarks
              if (abs (reader.hepeup.IDUP.at (iPart)) < 7) 
                {
                  v_f_quarks.push_back (dummy) ;        
                } // quarks
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 11)
                {
                  v_f_leptons_e.push_back (dummy) ;
                  idx_Ws_e.push_back(reader.hepeup.MOTHUP.at (iPart).first);
                }
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 13)
                {
                  v_f_leptons_mu.push_back (dummy) ;
                  idx_Ws_mu.push_back(reader.hepeup.MOTHUP.at (iPart).first);
                }
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 14 ||
                       abs (reader.hepeup.IDUP.at (iPart)) == 16)
                {
                  v_f_neutrinos.push_back (dummy) ;        
                }
              else if (abs (reader.hepeup.IDUP.at (iPart)) == 22 )
                {
                  v_f_photons.push_back (dummy) ;    
                }
            } // outgoing particles
        } // loop over particles in the event
      int warnNum = 0 ;
      if (v_f_quarks.size () < 0)
        {
          cout << "warning, not enough quarks" << endl ;
          ++warnNum ;
        }
      if (v_f_leptons_e.size ()+v_f_leptons_mu.size ()  >= 2)
        {
      //    cout << "warning, not enough leptons" << v_f_leptons_e.size ()+v_f_leptons_mu.size ()  << endl ;
          ++warnNum ;
        }
      if (v_f_leptons_mu.size () > 1)
        {
          //cout << "not e mu channel N_mu:" <<v_f_leptons_mu.size () << endl ;
          ++warnNum ;
        }
      if (v_f_neutrinos.size () > 2)
        {
          cout << "warning, to much neutrinos" << endl ;
          ++warnNum ;
        }
      if (warnNum > 0) continue ;

      
      //PG apply basic selections
      //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

     // cout << "!!! pass -----> N_mu: " <<v_f_leptons_mu.size ()<< " ---------->N_e: "<< v_f_leptons_mu.size () <<"------!!!" << endl ;
      //cout << "warning, not enough leptons" << v_f_leptons_e.size ()+v_f_leptons_mu.size ()  << endl ;
      TLorentzVector l1,l2,n1;
      if (v_f_leptons_e.size () >= 1)  l1=v_f_leptons_e.at(0);
      if (v_f_leptons_mu.size () >= 1) l2=v_f_leptons_mu.at(0);
      if (v_f_neutrinos.size () >= 1) n1=v_f_neutrinos.at(0);
      float weight = reader.hepeup.XWGTUP ;
      Histos.increaseNorm (weight) ;


      TLorentzVector ME = v_f_neutrinos.at (0)  ;
      
      float ptl1,ptl2 ;
      if (v_f_leptons_e.size () >= 1) ptl1= v_f_leptons_e.at (0).Pt () ;
      if (v_f_leptons_mu.size () >= 1) ptl2 = v_f_leptons_mu.at (0).Pt () ;
      if (ptl1 <10) continue;
      if (ptl2 <10) continue;

      Histos.fill ("no_shape", 0.5, weight) ;
      Histos.fill ("etal1", l1.Eta(), weight) ;
      Histos.fill ("etal2", l1.Eta(), weight) ;

      Histos.fill ("m_ptl1", ptl1, weight) ;
      Histos.fill ("m_ptl2", ptl2, weight) ;
      

      Histos.fill ("m_met", ME.Pt (), weight) ;
         

      if (max > 0 && max < events) 
        {
          cout << max << " events reached, exiting" << endl ;
          break ;
        }

    } //PG loop over input events

  return events ;
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


struct inputInfo
{
  double mass ;
  string mg_file ;
  double mg_xs ;
  string ph_file ;
  double ph_xs ;
  void print ()
    {
      cout << "----------------------------\n" ;
      cout << "mass    : " << mass << "\n" ;
      cout << "mg_file : " << mg_file << "\n" ;
      cout << "mg_xs   : " << mg_xs << "\n" ;
      cout << "ph_file : " << ph_file << "\n" ;
      cout << "ph_xs   : " << ph_xs << "\n" ;
      cout << "----------------------------\n" ;
    }
} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char ** argv) 
{

  if (argc < 2)
    {
      cerr << "Forgot to put the cfg file --> exit " << endl ;
      return 1 ;
    }

  int maxEventsPerSample = -1 ;
  if (argc >= 3)
    {
      maxEventsPerSample = atoi (argv[2]) ;
    }
  cout << "reading " << maxEventsPerSample << " events per sample" << endl ;

  CfgParser * gConfigParser = new CfgParser (argv[1]) ;

  //PG get the samples to be analised, 
  //PG including collection of LHE files and the relative XS

  vector<string> collections = gConfigParser->readStringListOpt ("general::samples") ;
  map<string, pair<float, vector<string> > > samples ;
  // loop over samples
  for (int i = 0 ; i < collections.size () ; ++i)
    {
      float XS = gConfigParser->readFloatOpt (collections.at (i) + "::XS") ;
      vector<string> inputfiles = gConfigParser->readStringListOpt (collections.at (i) + "::files") ;
      samples[collections.at (i)] = pair<float, vector<string> > (XS, inputfiles) ;
    } // loop over samples

  //PG prepare the histogram structs to be filled,
  //PG fill the histograms looping on LHE events

  // loop over samples
  map<string, pair<float, vector<string> > >::iterator it ;
  vector<histos> Histos ;
  for (it = samples.begin () ; it != samples.end () ; ++it)
    {
      cout << "sample: " << it->first << endl ;
      Histos.push_back (histos (it->first, it->second.first)) ;
      // loop over files
      int events = 0 ;
      for (int ifile = 0 ; ifile < it->second.second.size () ; ++ifile)
        {
          cout << "    reading: " << it->second.second.at (ifile) << endl ;
          std::ifstream ifs (it->second.second.at (ifile).c_str ()) ;
          LHEF::Reader reader (ifs) ;
          events += fillHistos (reader, Histos.back (), maxEventsPerSample) ;
          if (maxEventsPerSample > 0 && events >= maxEventsPerSample) break ;
          cout << "    read: " << events << " events" << endl ;

        } // loop over files
      Histos.back ().norm () ;
    } // loop over samples

  //PG save histograms
  string outfilename = gConfigParser->readStringOpt ("general::outputFile") ;
  TFile outFile (outfilename.c_str (), "recreate") ;
  for (int i = 0 ; i < Histos.size () ; ++i) Histos.at (i).save (outFile) ;
  outFile.Close () ;

  return 0 ;
}
