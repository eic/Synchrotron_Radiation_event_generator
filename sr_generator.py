import numpy as np
import pandas as pd
import yaml
import os
import sys
import ROOT
import matplotlib.pyplot as plt
import time
import argparse
#import pyhepmc_ng as hep
import pyhepmc as hep
import math

# =============================================================
class sr_generator:
    # ---------------------------------------
    def __init__(self,config):
        with open(config,'r') as stream:
            self.config = yaml.safe_load(stream)
          
        print('********************************************************')
        print('Synchrotron-radiation event generator from Synrad+ table')
        print('authors: Reynier Cruz-Torres*, Benjamen Sterwerf**')
        print('* Lawrence Berkeley National Laboratory')
        print('** University of California, Berkeley')
        print('')
        if 'input_single_photons' in self.config:
            self.input_single_photons = self.config['input_single_photons']
            print('Loading single photons from:',self.input_single_photons)
        else:
            print("Provide 'input_single_photons' in config file")
            exit()

        if 'n_events' in self.config:
            self.n_events = self.config['n_events']
            print("Number of events to run:",self.n_events)
        else:
            print("Number of events not provided in config file. Defaulting to 10")
            self.n_events = 10
            
        if 'integration_window' in self.config:
            self.integration_window = (float)(self.config['integration_window'])
            print("Time integration window selected:",self.integration_window)
        else:
            print("Time integration window not provided in config file. Defaulting to 100 ns")
            self.integration_window = (float)(1e-07)
            
        if 'seed' in self.config:
            self.seed = self.config['seed']
        else:
            self.seed = 0

        if self.seed == 0:
            print("The generator seed won't be constrained")
        else:
            print("The generator seed will be:",self.seed)
            ROOT.gRandom.SetSeed(self.seed)

        self.outpath = 'output_plots'
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        self.pid = 22 # 22 = photon, can be changed for testing
        self.status_000 = 4 # try 1, 2, or 4 (see page 31 in https://arxiv.org/pdf/1912.08005.pdf)
        self.status_xyz = 1

        if self.seed > 0:
            self.outfname = 'SR_out_int_window_{}ns_nevents_{}_pid_{}_status_{}_{}_seed_{}.hepmc'.format(self.integration_window*1e+09,
            self.n_events,self.pid,self.status_000,self.status_xyz,self.seed)
        else:
            self.outfname = 'SR_out_int_window_{}ns_nevents_{}_pid_{}_status_{}_{}.hepmc'.format(self.integration_window*1e+09,
            self.n_events,self.pid,self.status_000,self.status_xyz)

        if self.integration_window == 0:
            print("Time integration window is set to 0, writing out single photon events only.")
            print("IMPORTANT: We randomize WITHOUT replacing here! This is purely a hepmc conversion of a permutation of the original photons.")
            self.outfname = 'SR_out_single.hepmc'

        self.saveHisto = False
        if self.saveHisto :
            self.outHistFileName = self.outfname
            self.outHistFileName = self.outHistFileName.replace(".hepmc",".hist.root")

        print('Generated events will be saved in',self.outfname)
        self.f = hep.io.WriterAscii(self.outfname)

        print('********************************************************')

        self.load_single_photons()

    # ---------------------------------------
    def load_single_photons(self):
        print('Loading single photons')

        self.df = pd.read_csv(self.input_single_photons)
        self.df = self.df.sort_values('NormFact')
        # print(self.df[392:400])
        # sys.exit()

        n_entries = len(self.df)
        self.h1_df = ROOT.TH1D('h1_df',';entry;W [#gamma/sec]',n_entries,0,n_entries)
        hE = ROOT.TH1D('hE',';E [keV]',200,0,100)
        hT = ROOT.TH1D('hT',';t [ns]',200,0,100)
        hFlux = ROOT.TH1D('hFlux',';flux [#gamma/ns]',200,0,100)
        for i in range(n_entries):
            self.h1_df.SetBinContent(i+1,self.df['NormFact'].iloc[i])
            px = self.df['px'].iloc[i]
            py = self.df['py'].iloc[i]
            pz = self.df['pz'].iloc[i]
            E = math.sqrt( px**2 + py**2 + pz**2 )*1E6
            hE.Fill( E )
            nf = self.df['NormFact'].iloc[i] * 1e-9 # Hz -> GHz
            hFlux.Fill(nf)
            T = 1./nf if nf > 0 else 0
            hT.Fill( T )

        # plt.figure()
        # plt.hist(hE)
        # plt.savefig("hE.png",dpi=400)
        outHistFile = ROOT.TFile.Open ( "hE.root" ,"RECREATE")
        hE.Write()
        hT.Write()
        hFlux.Write()
        outHistFile.Close()
        sys.exit()

        if self.integration_window == 0 :
            if self.n_events > n_entries :
                print("Writing out a single photon permutation but you requested more single photons than arte available.")
                print("Setting n_events to {}".format(n_entries))
                self.n_events = n_entries

        # Prepare a pool of indices from which we can draw without replacement
        # We shuffle so that smaller samples retrieved with "mc head" aren't biased
        rng = np.random.default_rng(self.seed)
        self.lotto_bowl = list (range( 0,n_entries)) # see here for why a list https://stackoverflow.com/questions/20484195/typeerror-range-object-does-not-support-item-assignment
        rng.shuffle ( self.lotto_bowl )
         
    # ---------------------------------------
    def generate_an_event(self):
        event = []
        lastbin = self.h1_df.GetNbinsX()
        if self.integration_window==0 :
            x = self.lotto_bowl.pop() # get and remove the last entry
            photon = self.df.iloc[x]
            event.append(photon)
            return event,photon['NormFact']

        integrated_so_far = 0.
        while integrated_so_far < self.integration_window:
            x = self.h1_df.FindBin(self.h1_df.GetRandom())
            if x >= lastbin: continue
            photon = self.df.iloc[x]
            integrated_so_far += 1./photon['NormFact']
            event.append(photon)

        return event,1

    # ---------------------------------------
    def generate(self):
        print('Generating SR events')

        # For visualization
        events = []
        photons_per_event = []
        z_dist = []
        rho_dist = []

        ri = hep.GenRunInfo()
        ri.weight_names = ["NormFactor"]
        for i in range(self.n_events):
            if  i % 10000 == 0 : print('Working on event {} / {}'.format( i+1, self.n_events))

            event,NormFact = self.generate_an_event()
            # ---------------------------------------------------
            # Save to hepmc format
            # implemented following the example from:
            # https://github.com/scikit-hep/pyhepmc/blob/master/tests/test_basic.py
         
            evt = hep.GenEvent(hep.Units.GEV, hep.Units.MM)
            evt.run_info = ri
            if self.integration_window == 0: # could keep it for all events, this is just to keep the original output identical
                evt.weights=[NormFact]
            particles_out = []
            particles_in = []
            vertices = []

            # loop over each photon in the event
            for g in range(len(event)):
                x = event[g]['x']
                y = event[g]['y']
                z = event[g]['z']

                z_dist.append(z)
                rho_dist.append(np.sqrt(x*x+y*y))

                px = event[g]['px']
                py = event[g]['py']
                pz = event[g]['pz']
                E = np.sqrt(px**2 + py**2 + pz**2)

                pinx = E*x/np.sqrt(x*x+y*y+z*z)
                piny = E*y/np.sqrt(x*x+y*y+z*z)
                pinz = E*z/np.sqrt(x*x+y*y+z*z)

                # Particles going into the vertex
                pin = hep.GenParticle((pinx,piny,pinz,E),self.pid,self.status_000)
                pin.generated_mass = 0.
                evt.add_particle(pin)
                particles_in.append(pin)

                # Particles coming out of the vertex
                pout = hep.GenParticle((px,py,pz,E),self.pid,self.status_xyz)
                pout.generated_mass = 0.
                evt.add_particle(pout)
                particles_out.append(pout)

                # make sure vertex is not optimized away by the Writer
                v1 = hep.GenVertex((x,y,z,0.))
                v1.add_particle_in(pin)
                v1.add_particle_out(pout)
                evt.add_vertex(v1)
                vertices.append(v1)

            # ---------------------------------------------------
            self.f.write_event(evt)
            photons_per_event.append(len(event))

            # --------------------
            # for plotting
            if i < 6:
                events.append(event)

        # --------------------------
        # Make some plots
        self.plot_histo(photons_per_event,'# photons per event','Nphotons_per_event.png')
        self.plot_2d_scatter(events,'x','y','x [mm]','y [mm]',(-40,40),(-40,40),'events_x_v_y.png')
        self.plot_2d_scatter(events,'z','x','z [mm]','x [mm]',(-5000,3000),(-40,40),'events_z_v_x.png')
        self.plot_2d_scatter(events,'z','y','z [mm]','y [mm]',(-5000,3000),(-40,40),'events_z_v_y.png')

        # ----------------------------
        # Save histogram for later use
        if self.saveHisto :
            outHistFile = ROOT.TFile.Open ( self.outHistFileName ,"RECREATE")
            self.h1_df.Write()
            outHistFile.Close()

    # ---------------------------------------
    def plot_histo(self,hist,xlabel,fname):
        plt.figure()
        plt.hist(hist)
        plt.xlabel(xlabel)
        plt.ylabel('Counts')
        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath,fname),dpi=400)

    # ---------------------------------------
    def plot_2d_scatter(self,events,xl,yl,labelx,labely,xlim,ylim,fname):
        plt.figure(figsize=(13,8))
        for i in range(6):
            plt.subplot(2,3,i+1)
            x, y = [], []
            for j in range(len(events[i])):
                x.append(events[i][j][xl])
                y.append(events[i][j][yl])

        plt.scatter(x,y,marker='o',alpha=0.2)
        plt.xlim(xlim[0],xlim[1])
        plt.ylim(ylim[0],ylim[1])
        plt.xlabel(labelx)
        plt.ylabel(labely)
         
        plt.text(-8,0,r'{} $\gamma$'.format(len(events[i])),fontsize=15)

        plt.tight_layout()
        plt.savefig(os.path.join(self.outpath,fname),dpi=400)

# =============================================================
if __name__=='__main__':
    t0 = time.time()

    parser = argparse.ArgumentParser(description='Generating synchrotron radiation events')
    parser.add_argument('-c', '--configFile',
                            help='Path of config file for analysis',
                            action='store', type=str,
                            default='config.yaml')
    args = parser.parse_args()

    generator = sr_generator(args.configFile)
    generator.generate()

    print('Overall running time:',np.round((time.time()-t0)/60.,2),'min')
