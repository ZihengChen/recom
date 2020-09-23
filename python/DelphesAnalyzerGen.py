

from DelphesUltility import *

class DelphesAnalyzerGen():
    def __init__(self, fileName):
        self.events = uproot.open(fileName)["Delphes"]
        self.features = ['Event.Number']
        self.features+= ['Particle'+v for v in ['_size','.M1','.PID','.Status','.PT', '.Eta', '.Phi', '.Mass']]
        
    def run(self):
        
        out_gen = defaultdict(list) # out dict 

        # begin event loop
        for evts in self.events.iterate(self.features, namedecode="utf-8"):
            nev = len(evts['Event.Number'])
            for iev in range(nev):
                # process event iev

                ############################################################
                # object selection
                ############################################################
                
                # genpart loop
                genLeps, genPhotons, genGluons = [], [], []
                pids = [] # temp list for querying the mother pid
                for i in range(evts['Particle_size'][iev]): 
                    if len(genLeps) >= 2 and len(genPhotons) >= 1 and len(genGluons) >= 2: 
                        break 

                    pid  = evts['Particle.PID'][iev][i]
                    mIdx = evts['Particle.M1'][iev][i]
                    pids.append(pid)
                        
                    #find electron, anti-electron, photon and gluons
                    if (abs(pid)==11 or abs(pid)==13) and pids[mIdx] == 23: 
                        genLeps.append(i)
                    if pid==22 and (pids[mIdx]==3200211 or pids[mIdx]==3200112): # allow offshell G
                        genPhotons.append(i)
                    if pid==21 and (pids[mIdx]== 3200211 or pids[mIdx]==3160111): # allow offshell G
                        genGluons.append(i)
                    

                genLeps.sort(key=lambda i: evts['Particle.PT'][iev][i], reverse=True)
                genPhotons.sort(key=lambda i: evts['Particle.PT'][iev][i], reverse=True)
                genGluons.sort(key=lambda i: evts['Particle.PT'][iev][i], reverse=True)


                ############################################################
                # event selection
                ############################################################
                # should save all events
                if len(genLeps) >= 2 and len(genPhotons) >= 1 and len(genGluons) >= 2:
                    # no event cut
                    out_gen['eventNumber'].append(evts['Event.Number'][iev][0])
                    out_gen['nMuons'].append(len(genLeps)) # dummy 
                    out_gen['nElectrons'].append(len(genLeps)) # dummy
                    out_gen['nPhotons'].append(len(genPhotons)) 
                    out_gen['nJets'].append(len(genGluons))

                    # leading passing electron in the event
                    i,lep1 = genLeps[0], LorentzVector()
                    lep1.setptetaphim( evts['Particle.PT'][iev][i], evts['Particle.Eta'][iev][i], evts['Particle.Phi'][iev][i], evts['Particle.Mass'][iev][i] )
                    
                    # trailing passing electron in the event
                    i,lep2 = genLeps[1], LorentzVector()
                    lep2.setptetaphim( evts['Particle.PT'][iev][i], evts['Particle.Eta'][iev][i], evts['Particle.Phi'][iev][i], evts['Particle.Mass'][iev][i] )
                    
                    # leading passing photon in the event
                    i,gamma = genPhotons[0], LorentzVector()
                    gamma.setptetaphim( evts['Particle.PT'][iev][i], evts['Particle.Eta'][iev][i], evts['Particle.Phi'][iev][i], evts['Particle.Mass'][iev][i] )
                    
                    # leading passing jet in the event
                    i,jet1 = genGluons[0], LorentzVector()
                    jet1.setptetaphim( evts['Particle.PT'][iev][i], evts['Particle.Eta'][iev][i], evts['Particle.Phi'][iev][i], evts['Particle.Mass'][iev][i] )
                    
                    # Trailing passing jet in the event
                    i,jet2 = genGluons[1], LorentzVector()
                    jet2.setptetaphim( evts['Particle.PT'][iev][i], evts['Particle.Eta'][iev][i], evts['Particle.Phi'][iev][i], evts['Particle.Mass'][iev][i] )

                    fill_kinematics(out_gen, lep1, lep2, gamma, jet1, jet2)
                else:
                    print( evts['Event.Number'][iev][0] )
                    print(genLeps, genGluons, genPhotons)

        self.out_gen = pd.DataFrame(out_gen) # out dict
        