#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/EvtGen.h"
#include "Pythia8Plugins/Pythia8Rivet.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <chrono> // for time measurement


using namespace Pythia8;


inline double invariantMass(const Vec4 &p1, const Vec4 &p2, const Vec4 &p3) {
    Vec4 totalP = p1 + p2 + p3;
    return totalP.mCalc();
}


bool findpKpiold(const Event &event,int i){
    int dau1 = event[i].daughter1();
        int dau2 = event[i].daughter2();
        vector daulist = event[i].daughterList();
        bool foundProton = false, foundKaon = false, foundPion = false;
        if (dau1 > 0 && dau2 > 0) {
          Vec4 pProton, pKaon, pPion;

          
          //std::cout << "Tochterpartikel-Indizes: ";
           //   for (int idx : daulist) {
            //      std::cout << pythia.event[idx].id() << " ";
             // }
              //std::cout << std::endl;
          // Loop over daughters to find p, K-, and pi+
          for (int j = dau1; j <= dau2; ++j) {
            if (event[j].id() == 2212) { // Proton
              pProton = event[j].p();
              foundProton = true;
            } else if (event[j].id() == -321) { // K-
              pKaon = event[j].p();
              foundKaon = true;
            } else if (event[j].id() == 211) { // pi+
              pPion = event[j].p();
              foundPion = true;
            }
          }
          }
      if(foundProton && foundKaon && foundPion && daulist.size()==3) {return true;}
      else return false;
          }

// Find End Products independent of between decays



std::tuple<bool, Vec4, Vec4, Vec4, std::vector<int>> findpKpi(const Event &event,int i){
    Vec4 pProton, pKaon, pPion;
    std::vector<int> parList;
    int dau1 = event[i].daughter1();
    int dau2 = event[i].daughter2();
    //std::cout << "Tochterpartikel-Indizes: " << event[i].id() << " " << event[dau1].id() << " " << event[dau2].id() <<std::endl;
    //vector daulist = event[i].daughterList();
    int nProton = 0, nKaon = 0, nPion = 0, nOther=0;
    //std::cout << "Tochterpartikel-Indizes: " << event[dau1].id() << " " << event[dau2].id() <<std::endl;
    if (dau1 > 0 && dau2 > 0) {
      //std::cout << "Tochterpartikel-Indizes: ";
       //   for (int idx : daulist) {
        //      std::cout << pythia.event[idx].id() << " ";
         // }
          //std::cout << std::endl;
      // Loop over daughters to find p, K-, and pi+
      for (int j = dau1; j <= dau2; ++j) {
         parList.push_back(event[j].id());
         int dau3 = event[j].daughter1();
         int dau4 = event[j].daughter2();
             //std::cout << " and " << event[dau3].id() << " " << event[dau4].id() <<std::endl;
        if (dau3 > 0 && dau4 > 0){
            for (int k = dau3; k <= dau4; ++k) {
                  parList.push_back(event[k].id());
            if (event[k].id() == 2212) { // Proton
              pProton = event[k].p();
              nProton = nProton + 1;
            } else if (event[k].id() == -321) { // K-
              pKaon = event[k].p();
              nKaon = nKaon + 1;
            } else if (event[k].id() == 211) { // pi+
              pPion = event[k].p();
              nPion = nPion + 1;
            } else 
            {nOther = nOther + 1;}
            } 
            }
        else{
                    
        if (event[j].id() == 2212) { // Proton
          pProton = event[j].p();
          nProton = nProton + 1;
        } else if (event[j].id() == -321) { // K-
          pKaon = event[j].p();
          nKaon = nKaon + 1;
        } else if (event[j].id() == 211) { // pi+
          pPion = event[j].p();
          nPion = nPion + 1;
        } else 
            {nOther = nOther + 1;}
        } 
        }
      }
      
      // nOther==0
  if(nProton == 1 && nKaon ==1 && nPion == 1 && nOther==0) {return std::make_tuple(true,pProton,pKaon,pPion,parList);}
  else return std::make_tuple(false,pProton,pKaon,pPion,parList);
      }
 
int main(int argc, char* argv[]) {

  if (argc != 8) {
    cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide the arguments \n"
         << " 1. EvtGen decay file (e.g. DECAY_2010.DEC) \n"
         << " 2. EvtGen particle data (e.g. evt.pdl) \n"
         << " 3. PYTHIA8DATA path \n"
         << " 4. Flag to use EvtGen (true or false) \n"
         << " 5. Flag to use Rivet (true or false) \n"
	 << " 6. nEvents \n"
         << " 7. Output Name \n"
         << " Program stopped. " << endl;
    return 1;
  } 

  int nEvents = atoi(argv[6]); // Anzahl der Ereignisse
  int printInterval = 1000; // Intervall für die Fortschrittsanzeige
  string name = argv[7];
  string Rootend = ".root";
  string Plotend = ".pdf";

  //TString Rootbeg = strcat(Rootdir,name);
  TString Pythiaoutputname = name + Rootend;

  //TString Plotbeg = strcat(Plotdir,name);
  TString Plotname = name + Plotend;
  //TString Pythiaoutputname= "RootFiles/PythiaDatatest_new.root";
  //TString Plotname = "Plots/Plottest_new.pdf";
  // Check command-line arguments. 
  
  bool useEvtGen = (string(argv[4]) == "true");
  bool useRivet = (string(argv[5]) == "true");

  // Initialize Pythia.
  Pythia pythia; // col 13000
  pythia.readString("Beams:eCM = 13000.");          // Set collision energy
  pythia.readString("Beams:idA = 2212");            // Proton beam A
  pythia.readString("Beams:idB = 2212");            // Proton beam B
  pythia.readString("HardQCD:all = on");            // Enable hard QCD processes
  pythia.readString("PhaseSpace:pTHatMin = 20.");   // Minimum transverse momentum

  if (!useEvtGen) {
    cout << "Not using EvtGen." << endl;
    //pythia.readString("4122:onMode = off");         // Turn off all decays of Lambda_c+
    //pythia.readString("4122:onIfMatch = 2212 -321 211"); // Enable only Lambda_c+ -> p K- pi+
    pythia.readString("4122:mayDecay = on");        // Allow decays of Lambda_c+
    pythia.readString("4122:oneChannel = 1 1.0 0 2212 -321 211"); // Lambda_c+ -> p K- pi+
  } else {
    cout << "Using EvtGen." << endl;
  }

  // If Pythia fails to initialize, exit with error.
  if (!pythia.init()) return 1;

  // Initialize EvtGen if requested.
  EvtGenDecays *evtgen = nullptr;
  if (useEvtGen) {
    setenv("PYTHIA8DATA", argv[3], 1);
    evtgen = new EvtGenDecays(&pythia, argv[1], argv[2]);
    evtgen->readDecayFile(argv[1]);
  }

  // Create a Pythia8Rivet object and add (one or several) analyses.
  Pythia8Rivet rivet(pythia, name + ".yoda");
  if (useRivet) {
    rivet.addAnalysis("LHCB_2023_I2683025");
  }

  
  // ROOT-Histogramme erstellen
    TH1D *hInvMass = new TH1D("hInvMass", "Invariant Mass of #Lambda_{c}^{+} -> p K^{-} #pi^{+};Mass [GeV];Entries", 10000, 2.28645999999, 2.28646000001);
    TH1D *hpkMass = new TH1D("hpkMass", "Invariant Mass of p K^{-};Mass [GeV];Entries", 10000, 1.0, 2.5);
    TH1D *hkpiMass = new TH1D("hkpiMass", "Invariant Mass of K^{-} #pi^{+};Mass [GeV];Entries", 10000, 0., 1.5);
    TH1D *hkMass = new TH1D("hkMass", "Invariant Mass of K^{-};Mass [GeV];Entries", 10000, 0., 1.5);
    TH1D *hppiMass = new TH1D("hppiMass", "Invariant Mass of p #pi^{+};Mass [GeV];Entries", 10000, 0., 2.5);
    TH1D *hpiMass = new TH1D("hpiMass", "Invariant Mass of #pi^{+};Mass [GeV];Entries", 10000, 0., 1.5);
    TH1D *hpMass = new TH1D("hpMass", "Invariant Mass of p^{+};Mass [GeV];Entries", 10000, 0., 1.5);
    TH2D *hDalitz = new TH2D("hDalitz", "Dalitz Plot;M^{2}(pK^{-}) [GeV^{2}];M^{2}(K^{-}#pi^{+}) [GeV^{2}]", 200, 1.5, 5.0, 200, 0.0, 2.0);

    // ROOT-Datei zum Speichern der Histogramme öffnen
    TFile *outputFile = new TFile(Pythiaoutputname, "RECREATE");
    if (!outputFile->IsOpen()) {
        std::cerr << "Error: Unable to create ROOT file!" << std::endl;
        return 1;
    }
    
  // Create histograms for invariant masses.
  //Hist invariantMassHist("Invariant Mass of Lambda_c+ -> p K- pi+", 100, 2.2, 2.4);
  //Hist pkMassHist("Invariant Mass of pK+",100,1.5,2.2);

  

  // Startzeit erfassen
  auto startTime = std::chrono::steady_clock::now();

  // Event loop.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;

    // Perform decays with EvtGen if enabled.
    if (evtgen) evtgen->decay();
 // Fortschrittsanzeige alle 'printInterval' Ereignisse
        if (iEvent % printInterval == 0 && iEvent > 0) {
            // Aktuelle Zeit erfassen
            auto currentTime = std::chrono::steady_clock::now();
            double elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double estimatedTotalTime = (elapsedTime / iEvent) * nEvents;
            int remainingTime = estimatedTotalTime - elapsedTime;
                        int remainhours = remainingTime/3600;
			int remainmin = remainingTime/60 - remainhours*60;
			int remainsec = remainingTime - remainhours*3600- remainmin*60;

            std::cout << "Generierte Ereignisse: " << iEvent
                      << " / " << nEvents
                      << " (" << (100.0 * iEvent / nEvents) << "% abgeschlossen), "
                    //  << "verbleibende Zeit: ~" << remainingTime << " Sekunden" << "or" << remainmin << ":" << remainsec
	                  << "verbleibende Zeit: ~" << remainhours << ":" << remainmin << ":" << remainsec	  
                      << std::endl;
        }

    
    // Loop over particles in the event.
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 4122) { // Find Lambda_c+
       /*
        int dau1 = pythia.event[i].daughter1();
        int dau2 = pythia.event[i].daughter2();

       
        vector daulist = pythia.event[i].daughterList();
        if (dau1 > 0 && dau2 > 0) {
          Vec4 pProton1, pKaon1, pPion1;
          bool foundProton = false, foundKaon = false, foundPion = false;
          
          //std::cout << "Tochterpartikel-Indizes: ";
           //   for (int idx : daulist) {
            //      std::cout << pythia.event[idx].id() << " ";
             // }
              //std::cout << std::endl;
          // Loop over daughters to find p, K-, and pi+
          for (int j = dau1; j <= dau2; ++j) {
            if (pythia.event[j].id() == 2212) { // Proton
              pProton = pythia.event[j].p();
              foundProton = true;
            } else if (pythia.event[j].id() == -321) { // K-
              pKaon = pythia.event[j].p();
              foundKaon = true;
            } else if (pythia.event[j].id() == 211) { // pi+
              pPion = pythia.event[j].p();
              foundPion = true;
            }
          }
          }
          if (foundProton && foundKaon && foundPion && daulist.size()==3 ) {
          */
          auto [success, pProton, pKaon, pPion,daulist] = findpKpi(pythia.event,i);
          

          // Wenn alle Partikel gefunden wurden, berechne die Massen
          if (success) {
              double invMass = invariantMass(pProton, pKaon, pPion);
              hInvMass->Fill(invMass);
              
              double pkMass = invariantMass(pProton, pKaon, 0);
              double kpiMass = invariantMass(pKaon, pPion, 0);
              double ppiMass = invariantMass(pProton, pPion, 0);
              hpkMass->Fill(pkMass);
              hkpiMass->Fill(kpiMass);
              hppiMass->Fill(ppiMass);
              hkMass->Fill(pKaon.mCalc());
              hpiMass->Fill(pPion.mCalc());
              hpMass->Fill(pProton.mCalc());
              // Berechnung der invarianten Massenquadrate für den Dalitz-Plot
              double m2_pK = pow(invariantMass(pProton, pKaon, 0), 2);
              double m2_Kpi = pow(invariantMass(pKaon, pPion, 0), 2);
              hDalitz->Fill(m2_pK, m2_Kpi);
              
              //std::cout << "Inhalt der Liste: ";
              //for (int value : daulist) {
              //    std::cout << value << " ";
              //    }
              //std::cout << std::endl;
          /*
              // Tochterliste ausgeben
              std::cout << "Tochterpartikel-Indizes: ";
              for (int idx : daulist) {
                  std::cout << pythia.event[idx].id() << " ";
              }
              std::cout << std::endl;*/
          
        }
      }
    }
    // Push event to Rivet.
    if (useRivet) rivet();
  }

  
  // Statistiken ausgeben
  pythia.stat();

  // Tell Rivet to finalise the run.
  if (useRivet) rivet.done();

  // ROOT-CANVAS für die Ausgabe
  TCanvas *c1 = new TCanvas("c1", "Invariant Mass and Dalitz Plot", 1200, 600);
  c1->Divide(2, 2); // Zwei nebeneinander liegende Plots

  // Histogramm für die invariante Masse zeichnen
  c1->cd(1);
  hInvMass->SetLineColor(kBlue);
  hInvMass->Draw();

  // Dalitz-Plot zeichnen
  c1->cd(2);
  gStyle->SetOptStat(0);
  hDalitz->Draw("COLZ");
  
  c1->cd(3);
  hpkMass->SetLineColor(kBlue);
  hpkMass->Draw();

  c1->cd(4);
  hkpiMass->SetLineColor(kBlue);
  hkpiMass->Draw();
  
  // Canvas speichern und anzeigen
  c1->SaveAs(Plotname);

  // Histogramme in der ROOT-Datei speichern
  outputFile->cd();
  hInvMass->Write();
  hDalitz->Write();
  hpkMass->Write();
  hkpiMass->Write();
  hppiMass->Write();
  hpMass->Write();
  hkMass->Write();
  hpiMass->Write();
  

  // ROOT-Datei schließen
  outputFile->Close();
  std::cout << "Daten wurden in "<< Pythiaoutputname << " gespeichert." << std::endl;

  // Speicher freigeben
  delete hInvMass;
  delete hDalitz;
  delete outputFile;
  return 0;
}

