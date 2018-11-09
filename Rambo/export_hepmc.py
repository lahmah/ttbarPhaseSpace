import numpy as np
import os

def export_hepmc(E_CM, data, weights, filename):
    n_out = int(data.shape[1]/4)
    
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w") as f:
        f.write("\nHepMC::Version 2.06.09\nHepMC::IO_GenEvent-START_EVENT_LISTING\n")
        for i in range(data.shape[0]):
            # event
            # E evt_number no_mpi scale alphq_qcd alpha_qed signal_id barcode_signal_process_vertex no_vertices barcode_particle_1 barcode_particle_2 no_random_state {random_state} no_weights {weights}
            f.write("E %i -1 0 1.0000000000000000e+00 1.0000000000000000e+00 0 0 1 10001 10002 0 1 %e\n" % (i, weights[i]))
            
            # weights
            f.write("N 1 \"0\"\n")
            
            # units
            f.write("U GEV MM\n")
            
            # vertex
            # V barcode id x y z ctau no_incoming no_outgoing no_weights {weights}
            f.write("V -1 0 0 0 0 0 2 %i 0\n" % n_out)
            
            # incoming particles
            # P barcode PDG_id px py pz energy gen_mass status_code pol_theta pol_phi barcode_vertex_incoming no_flow {code_index, code}
            f.write("P 10001 11 0 0 %e %e 0 4 0 0 -1 0\n" % (E_CM/2, E_CM/2))
            f.write("P 10002 -11 0 0 %e %e 0 4 0 0 -1 0\n" % (-E_CM/2, E_CM/2))
            
            # outgoing particles
# 6 -> 5 24 ->       5 -13 14
# -6 -> -5 -24 ->    -5 11 -12
            for j in range(n_out):
                if j == 0:
                    pid = -24 
                elif j == 1:
                    pid = 6
                elif j == 2:
                    pid = -5 
               
                E = data[i, 4*j]
                px = data[i, 4*j+1]
                py = data[i, 4*j+2]
                pz = data[i, 4*j+3]
                f.write("P %i %i %e %e %e %e 0 1 0 0 0 0\n" % (10003+j, pid, px, py, pz, E))
            
        f.write("HepMC::IO_GenEvent-END_EVENT_LISTING")
