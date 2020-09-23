

from DelphesUltility import *



class DelphesAnalyzer():
    def __init__(self, fileName):
        self.events = uproot.open(fileName)["Delphes"]
        self.features = ['Event.Number','Photon*','Electron*']
        self.features+= ['Muon'+v for v in ['_size','.PT','.Eta','.Phi']]
        self.features+= ['Jet'+v for v in ['_size','.PT','.Eta','.Phi','.Mass',".Flavor"]]
        
    def run(self):
        
        out_eeg = defaultdict(list) # out dict 
        out_mmg = defaultdict(list) # out dict 

        # begin event loop
        for evts in self.events.iterate(self.features, namedecode="utf-8"):
            nev = len(evts['Event.Number'])
            for iev in range(nev):
                # process event iev

                ############################################################
                # object selection
                ############################################################
                
                # muon loop select good muons
                pass_muons = []
                for i in range(evts["Muon_size"][iev]):
                    if evts['Muon.PT'][iev][i] > 10 and abs(evts['Muon.Eta'][iev][i] < 2.4):
                        pass_muons.append(i)
                pass_muons.sort(key=lambda i: evts['Muon.PT'][iev][i], reverse=True)
                
                
                # electron loop select good electrons
                pass_electrons = []
                for i in range(evts["Electron_size"][iev]):
                    if evts['Electron.PT'][iev][i] > 20 and abs(evts['Electron.Eta'][iev][i] < 2.5):
                            pass_electrons.append(i)
                pass_electrons.sort(key=lambda i: evts['Electron.PT'][iev][i], reverse=True)
            
                # gamma loop: select good photon
                pass_photons = []
                for i in range(evts["Photon_size"][iev]):
                    if evts['Photon.PT'][iev][i] > 15 and abs(evts['Photon.Eta'][iev][i] < 2.3):
                        pass_photons.append(i)
                pass_photons.sort(key=lambda i: evts['Photon.PT'][iev][i], reverse=True)
                    
                # jet loop: select good jets
                pass_jets = []
                for i in range(evts["Jet_size"][iev]):
                    if evts['Jet.PT'][iev][i] > 30 and abs(evts['Jet.Eta'][iev][i] < 2.3):
                        pass_jets.append(i)
                pass_jets.sort(key=lambda i: evts['Jet.PT'][iev][i], reverse=True)
            
            


                ############################################################
                # event selection
                ############################################################

                # eeg channel
                if (len(pass_electrons)>=2 and len(pass_photons)>=1 and len(pass_jets)>=2):
                    # channel selection. m_ll in z window, leading pt >30 ..
                    leadinglepton_pt = evts['Electron.PT'][iev][pass_electrons[0]]
                    if leadinglepton_pt < 30:
                        continue
                
                    # fill event general information
                    out_eeg['eventNumber'].append(evts['Event.Number'][iev][0])
                    out_eeg['nMuons'].append(len(pass_muons))
                    out_eeg['nElectrons'].append(len(pass_electrons))
                    out_eeg['nPhotons'].append(len(pass_photons))
                    out_eeg['nJets'].append(len(pass_jets))
                    
                    # leading passing electron in the event
                    i,lep1 = pass_electrons[0], LorentzVector()
                    lep1.setptetaphim( evts['Electron.PT'][iev][i], evts['Electron.Eta'][iev][i], evts['Electron.Phi'][iev][i], 0.00051)
                    # trailing passing electron in the event
                    i,lep2 = pass_electrons[1], LorentzVector()
                    lep2.setptetaphim( evts['Electron.PT'][iev][i], evts['Electron.Eta'][iev][i], evts['Electron.Phi'][iev][i], 0.00051)
                    # leading passing photon in the event
                    i,gamma = pass_photons[0], LorentzVector()
                    gamma.setptetaphim(evts['Photon.PT'][iev][i], evts['Photon.Eta'][iev][i],evts['Photon.Phi'][iev][i], 0.0)
                    # leading passing jet in the event
                    i,jet1 = pass_jets[0], LorentzVector()
                    jet1.setptetaphim(evts['Jet.PT'][iev][i], evts['Jet.Eta'][iev][i], evts['Jet.Phi'][iev][i], evts['Jet.Mass'][iev][i])
                    # Trailing passing jet in the event
                    i,jet2 = pass_jets[1], LorentzVector()
                    jet2.setptetaphim(evts['Jet.PT'][iev][i], evts['Jet.Eta'][iev][i], evts['Jet.Phi'][iev][i], evts['Jet.Mass'][iev][i])

                    fill_kinematics(out_eeg, lep1, lep2, gamma, jet1, jet2)



                # mmg channel
                elif (len(pass_muons)>=2 and len(pass_photons)>=1 and len(pass_jets)>=2):

                    # channel selection. m_ee in z window, leading pt >30 ..
                    leadinglepton_pt = evts['Muon.PT'][iev][pass_muons[0]]
                    if leadinglepton_pt < 30:
                        continue
                
                    # fill event general information
                    out_mmg['eventNumber'].append(evts['Event.Number'][iev][0])
                    out_mmg['nMuons'].append(len(pass_muons))
                    out_mmg['nElectrons'].append(len(pass_electrons))
                    out_mmg['nPhotons'].append(len(pass_photons))
                    out_mmg['nJets'].append(len(pass_jets))
                    
                    # leading passing electron in the event
                    i,lep1 = pass_muons[0], LorentzVector()
                    lep1.setptetaphim( evts['Muon.PT'][iev][i], evts['Muon.Eta'][iev][i], evts['Muon.Phi'][iev][i], 0.1056)
                    # trailing passing electron in the event
                    i,lep2 = pass_muons[1], LorentzVector()
                    lep2.setptetaphim( evts['Muon.PT'][iev][i], evts['Muon.Eta'][iev][i], evts['Muon.Phi'][iev][i], 0.1056)
                    # leading passing photon in the event
                    i,gamma = pass_photons[0], LorentzVector()
                    gamma.setptetaphim(evts['Photon.PT'][iev][i], evts['Photon.Eta'][iev][i],evts['Photon.Phi'][iev][i], 0.0)
                    # leading passing jet in the event
                    i,jet1 = pass_jets[0], LorentzVector()
                    jet1.setptetaphim(evts['Jet.PT'][iev][i], evts['Jet.Eta'][iev][i], evts['Jet.Phi'][iev][i], evts['Jet.Mass'][iev][i])
                    # Trailing passing jet in the event
                    i,jet2 = pass_jets[1], LorentzVector()
                    jet2.setptetaphim(evts['Jet.PT'][iev][i], evts['Jet.Eta'][iev][i], evts['Jet.Phi'][iev][i], evts['Jet.Mass'][iev][i])

                    fill_kinematics(out_mmg, lep1, lep2, gamma, jet1, jet2)
        
        self.out_eeg = pd.DataFrame(out_eeg) # out dict 
        self.out_mmg = pd.DataFrame(out_mmg) # out dict 