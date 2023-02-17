import os
import numpy as np
import pandas as pd

def run_command(cmd):
    import subprocess
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output.decode('utf-8')
                
class RNAstructure(object): 
    def __init__(self, rnastructure_path) -> None:
        self.rnastructure_path = rnastructure_path if rnastructure_path[-1] == '/' else rnastructure_path+'/'
        self.directory = 'temp/rnastructure/'

    def fit(self, sequence, reference='reference'):
        self.sequence = sequence
        self.__make_temp_folder()
        self.__make_files()
        self.__create_fasta_file(reference, sequence)

    def predict_ensemble_energy(self):
        cmd = f"{self.rnastructure_path}EnsembleEnergy {self.fasta_file} --sequence"
        splitted_output = run_command(cmd).split(' ')
        return float(splitted_output[splitted_output.index(f"kcal/mol\n\nEnsemble")-1])

    def predict_partition(self, temperature_k =None):
        cmd = f"{self.rnastructure_path}partition {self.fasta_file} {self.pfs_file}"
        if temperature_k != None:
            cmd += ' --temperature '+str(temperature_k)
        run_command(cmd)
        run_command(self.rnastructure_path+'ProbabilityPlot '+ self.pfs_file + ' -t '+self.prob_file)
        with open(self.prob_file,"r") as f:
            lines=f.readlines()
            out={'i':[],'j':[],'p':[]}
            for x in range(len(lines)):
                if x>1:
                    ls=lines[x].split("\t")
                    out["i"]+=[int(ls[0])]
                    out["j"]+=[int(ls[1])]
                    out["p"]+=[float(ls[2])]
        return self.__cast_pairing_prob(out)

    def predict_reference_deltaG(self):
        cmd = f"{self.rnastructure_path}Fold {self.fasta_file} {self.ct_file}"
        run_command(cmd)
        assert os.path.getsize(self.ct_file) != 0, f"{self.ct_file} is empty, check that RNAstructure works"
        return self.__extract_deltaG_struct()

    def __make_temp_folder(self):
        isExist = os.path.exists(self.directory)
        if not isExist:
            os.makedirs(self.directory)
        return self.directory

    def __make_files(self, temp_prefix='temp'):
        self.pfs_file = f"{self.directory}/{temp_prefix}.pfs"
        self.ct_file = f"{self.directory}/{temp_prefix}.ct"
        self.dot_file = f"{self.directory}/{temp_prefix}_dot.txt"
        self.fasta_file = self.directory+'/'+temp_prefix+'.fasta'
        self.prob_file = self.directory+'/'+temp_prefix+'_prob.txt'

    def __create_fasta_file(self, reference, sequence):
        # push the ref into a temp file
        temp_fasta = open(self.fasta_file, 'w')
        temp_fasta.write('>'+reference+'\n'+sequence)
        temp_fasta.close()

    # cast the temp file into a dot_bracket structure and extract the attributes
    def __extract_deltaG_struct(self):
        run_command(f"ct2dot {self.ct_file} 1 {self.dot_file}")
        temp_dot = open(self.dot_file, 'r')
        first_line = temp_dot.readline().split()
        # If only dots in the structure, no deltaG 
        print('first_line',first_line)
        if len(first_line) == 4:
            _, _, deltaG, _ = first_line
            deltaG = float(deltaG)
        if len(first_line) == 1:
            deltaG, _ = 0.0, first_line[0][1:]

        sequence = temp_dot.readline()[:-1] #  Remove the \n
        structure = temp_dot.readline()[:-1] # Remove the \n
        return deltaG, structure

    def run(self, sequence):
        if os.path.exists(os.path.join(self.directory, sequence+'.txt')):
            # read from file
            with open(os.path.join(self.directory, sequence+'.txt'), 'r') as f:
                deltaG, structure = f.readline().split()
        else:
            self.fit(sequence)
            deltaG, structure = self.predict_reference_deltaG()
            with open(os.path.join(self.directory, sequence+'.txt'), 'w') as f:
                f.write(f"{deltaG} {structure}")
        return deltaG, structure

if __name__ == "__main__":
    rna = RNAstructure('/Users/ymdt/src/RNAstructure/exe')
    rna.fit(sequence='AAGATATTCGAAAGAATATCTT')
    print("DeltaG + structure:", rna.predict_reference_deltaG())
    print("Ens. energy:", rna.predict_ensemble_energy())
    print("One line command:", rna.run('AAGATATTCGAAAGAATATCTT'))