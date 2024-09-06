import re
import os, shutil
import subprocess
import glob
import time
import pandas as pd
import numpy as np
import biobox as bb
import logging

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.ResidueDepth import min_dist, get_surface, residue_depth
except:
    print("biopython and msms unavailable. Will not be able to calculate residue depth")


class Measure(object):
    
    def __init__(self, df_input, outdir="result", activate_log=False, log_path='measure_log.txt',
                 features=["propka", "pkaANI", "sasa", "depth"]):
        
        self.activate_log = False
        if activate_log:
            self.activate_log = activate_log
            self.log_path = os.path.join(outdir, log_path)

            # define a logger
            self.logger = logging.getLogger('MeasureLog')
            self.logger.setLevel(level = logging.DEBUG)

            formatter = logging.Formatter('%(message)s') # as simple as possible
            handler = logging.FileHandler(self.log_path, encoding = 'UTF-8')
            handler.setLevel(logging.INFO)
            handler.setFormatter(formatter)

            self.logger.addHandler(handler)

        self._setup_measures(features)

        # for restarting
        self.current_index = 0
        
        # document failed pdb files
        self.wrong_pdb_file = list()
        
        self.outdir = outdir
        self.df_input = df_input
        self.folder = os.path.join(outdir, "curated")
               
        # Check that all files in DataFrame appear at least once in folder
        files1=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))] #to find AlphaFold entries
        files2=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))] #to find PDB entries
        for f in df_input["PDB Code"].values:
            if f not in files1 and f not in files2:
                print("WARNING: %s not found in folder %s"%(f, self.folder))

        self.pkaoutdir = os.path.join(outdir, "propkaoutput")
        if not os.path.exists(self.pkaoutdir):
            os.makedirs(self.pkaoutdir)
        
        self.PDB_only = False
        
        if 'Uniprot_Entry' in self.df_input.columns: 
            columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid']
            self.df = pd.DataFrame(columns=columns)
        else:
            self.PDB_only = True
            columns = ['PDB Code', 'Chain', 'Resid', 'pKa', 'sasa', 'depth']
            self.df = pd.DataFrame(columns = columns)
                
    def _setup_measures(self, features):
        '''
        convert a list of features into a measuring protocol
        '''

        # measures to carry out [label for DataFrame column, and function evaluating a file]
        # functions must return a dataframe [chain, resid, measure]
        self.measures = []
        for m in features:
            
            if m == "propka":
                self.measures.append([m, self.calculate_pka])
            elif m == "pkaANI":
                self.measures.append([m, self.calculate_pkaANI])
            elif m == "sasa":
                self.measures.append([m, self.calculate_sasa])
            elif m == "depth":    
                self.measures.append([m, self.calculate_depth])
            else:
                raise Exception(f"measure {m} unknown")
    
    def save_state(self, outname="measures.csv"):
        '''
        Save a csv file in output directory
        '''
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)

    
    def measure_dataframe(self):
        
        # use a different method if handling pdb codes only
        if self.PDB_only:
            return 'Call PDB_only method'
        
        files = glob.glob(os.path.join(self.folder, "*pdb"))
        
        first_index = self.df_input.index[0]
        
        for index, row in self.df_input.iterrows():
            
           
            PDBCODE = row["PDB Code"]
            chains = row["Chains"].split("/")
            self.current_index = index
            
            uniprot_code = row["Uniprot_Entry"]
            print("\n# UNIPROT: %s PDB: %s, chain(s): %s"%(uniprot_code, PDBCODE, " ".join(chains)))
            
            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files:
                
                if PDBCODE not in f:
                    continue

                tstart = time.time()
                print("\n> File: %s"%f)
                             
                # create temporary DataFrame for data of current file,
                # to be then appended to main DataFrame self.df
                
                columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid']
                df = pd.DataFrame(columns=columns)
                
                
                # append to temporary DataFrame all lysines in the file of interest             
                try:
                    M = bb.Molecule(f) # sometimes bb does not work with a pdb file
                except:
                    self.wrong_pdb_file.append(f)
                    continue
                
                _, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
                for i in idxs:

                    #save only lysine entries from chain of interest
                    if M.data["chain"].values[i] not in chains:
                        continue
                    
                    
                    data = ({'Uniprot_Entry': uniprot_code,
                        'PDB Code': f.split(".")[0],
                        'Method': row["Method"],
                        'Resolution': row["Resolution"],
                        'Chain': M.data["chain"].values[i],
                        'Resid': M.data["resid"].values[i]})

                    df = pd.concat([df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                
                print(f">> {len(df)} lysines of interest found")

                # iterate over measures to carry out (according to self.measures)
                for meas in self.measures:
                    print(f">> evaluating {meas[0]}...")
                    try:
                        df[meas[0]] = np.nan # create new column for measure
                        result = meas[1](f) # run measurement
                        #print("result is ", result)
                        try:
                            df = self._combine_dataframes(df, result, meas[0]) #insert measures into temporary DataFrame
                        except Exception as e:
                            print("combination failed ", str(e))
                                
                    except Exception as e:
                        print(f"ERROR in evaluation: {e}")
                        continue
                
                print(">> file processed in %4.2f sec."%(time.time()-tstart))
                
                # document the data to a log file
                if self.activate_log:
                    if df.empty == False:
                        try:
                            self.logger.info(df.to_string().strip('    Uniprot_Entry                    PDB Code Method Resolution Chain Resid    pKa       sasa'))
                            self.logger.info('--------------------------------------------------------------------------')
                        except:
                            print('Error in logging.')
                
                #append temporary DataFrame with all measures on a single file to main DataFrame
                self.df = pd.concat([self.df, df], ignore_index=True)
                
            # remove possible duplicated rows (if restarted)
            if index == first_index:
                self.df = self.df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)

    
    def recover_from_log(self, log_path):
    
        if self.PDB_only:
            return 'Function not callable.'
        
        columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid', 'pKa', 'sasa']
        log_to_df = pd.DataFrame(columns=columns)
        
        with open(log_path) as inf:
            for line in inf:
                line = line.replace('--------------------------------------------------------------------------',' ')
                parts = line.split()
                if len(parts) == 0:
                    continue
                parts = parts[1:]
                data = dict(zip(columns, parts))
                log_to_df = log_to_df.append(data, ignore_index=True)
        
        return log_to_df
    
    
    def restart_measure(self):
        
        if self.PDB_only:
            return 'Function not callable'
            
        self.df_input = self.df_input.loc[self.current_index:, :]
        self.measure_dataframe()
      
  
    def _combine_dataframes(self, target, to_merge, col_name):
        '''
        target is a DataFrame to be filled with data, to_merge contains the data.
        Values to insert are indexed in both array by two columns: Chain and Resid.
        '''
        # e.g. self._combine_dataframes(df, result, meas[0])
        
        if col_name == "propka":
            col_name = "pKa"
        
        for i in range(len(target)):
            
            chain_value = target.loc[i, "Chain"]
            resid_value = target.loc[i, "Resid"]
            
            idx = np.where((to_merge["Chain"] == chain_value) & (to_merge["Resid"].astype(int) == resid_value))
            if len(idx[0]) == 0:
                continue
            
            #print("<<<<<<<<<", col_name)
            #print( to_merge)
            target.at[i, col_name] = to_merge.loc[idx[0][0], col_name]
            #print(">>>>>>>>>>")
          
        return target


    def measure_PDB_only(self):
        
        if not self.PDB_only:
            return 'Calling the wrong method.'
        
        files = glob.glob(os.path.join(self.folder, "*pdb"))
        
        for _, row in self.df_input.iterrows():
            PDBCODE = row['PDB Code']
            
            for f in files:
                if PDBCODE not in f:
                    continue
                    
                tstart = time.time()
                print("\n> File: %s"%f)
                
                result_list = list()
                for meas in self.measures:
                    print(">> evaluating %s..."%meas[0])
                    try:
                       
                        result = meas[1](f) # run measurement
                        result_list.append(result)
                        
                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue
                try:
                    # add sasa data to pka data
                    to_merge = result_list[1]
                    target = result_list[0]
                    
                    for i in range(len(to_merge)):      
    
                        chain_value = to_merge.loc[i, "Chain"]
                        resid_value = to_merge.loc[i, "Resid"]
    
                        idx = np.where((target["Chain"] == chain_value) & (target["Resid"] == resid_value))
    
                        if len(idx[0]) == 0:
                            row = to_merge.loc[i,:]
                            target = target.append(row, ignore_index = True)
                            continue
    
                        target.at[idx[0][0], 'sasa'] = to_merge.loc[i, 'sasa']
                   
                    # final thing to do: make sure the order of the columns is correct
                    target = target[['Chain', 'Resid', 'pKa', 'sasa']]
                    print(">> %s lysines of interest found"%len(target))
                    
                    # todo1: add a col i.e. PDB Code
                    target.insert(0, 'PDB Code', [PDBCODE] * len(target))
                    
                    # todo2: append to the main df
                    self.df = pd.concat([self.df, target], ignore_index=True)
                    
                except Exception as e:
                    print("ERROR in data aggregation: %s"%e)
                print(">> file processed in %4.2f sec."%(time.time()-tstart))             
    
   
    def calculate_pka(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data
        parse errors, and return a dataframe containing all measurements not yielding an error.
        '''
        
        code_for_df = os.path.basename(path).split(".")[0]
        error_file_name = os.path.join(self.pkaoutdir, "%s_propka_errors.txt"%code_for_df)
        #print('start propka calc')
        try:   
            f = open(error_file_name, 'w')
            process = subprocess.Popen(['python', '-m', 'propka', path],
                                stdout=f, stderr=f)
            stdout, stderr = process.communicate()
            f.close()
                    
        except Exception as e:
            print("propka call failed: ", str(e))
            f.close()
            
            try:
                shutil.move(code_for_df, os.path.join(self.pkaoutdir, code_for_df))
            except:
                pass
            
            raise Exception('Failed to obtain pKa data. %s.'%e)
 
        #print('begin parsing')   
 
        try:
            propka_lys_fails = self.parse_propka_errors(error_file_name)
        except Exception as e:
            raise Exception("Failed extracting PROPKA errors. %s"%e)
           
        try:
            pkafile = code_for_df + '.pka'
            propres = open(pkafile)
        except Exception:
            raise Exception('Failed to find %s'%pkafile)
                
        lys_number = list()
        pkas = list()
        chain = list() 
        try:
            for line in propres:
                if re.search('^   LYS' , line):
                    try:
    
                        line = line[6:]
                        line = line.split()

                        # reject adding entries associated with errors in structure (as per logfile)
                        if len(propka_lys_fails)>0:
                            idx = np.where((propka_lys_fails["Chain"] == line[1]) & (propka_lys_fails["Resid"].astype(int) == int(line[0])))
                            if len(idx[0])>0:
                                continue    
            
                        lys_number.append(int(line[0]))
                        chain.append(line[1])
                        pkas.append(float(line[2]))
                        
                    except Exception as e:
                        print("> Error %s"%e)
                        continue
    
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
        
        except Exception as e:
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
            raise Exception('Failure parsing %s.pka. %s'%(code_for_df, e))
            
        #print('assemblying dataframe')
        try:
            df = pd.DataFrame({'Resid':lys_number,
                               'Chain': chain,
                               'pKa':pkas})
            
            df.sort_values(by=['pKa'], inplace=True)
            df = df.dropna()
            df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
            
        except Exception as e:
            raise Exception('Failed to construct pKa dataframe. %s'%e)
    
        #print(df)
        return df
    

    def parse_propka_errors(self, path):
        '''
        parse the PROPKA output file and appends unique chain and resid of any lysines mentioned a DataFrame.
        This list is returned to main and later the residues in it are removed from the df.
        '''
        
        f = open(path, 'r')
        list_remove = list()
        cnt = 0
        for line in f:
            cnt += 1
            lys_raw = re.findall('LYS [\d]*[\s][\w]*', line)
    
            for line in lys_raw:
                words = line.split(' ')
                resid = words[1]
                chain = words[2]
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])  
    
            lys_raw_2 = re.findall('[\d]*-LYS \(\w\)', line)
    
            for line in lys_raw_2:
                words = line.split()
                chain = (words[1])[1:-1]
                words_2 = line.split('-')
                resid = words_2[0]
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])  
        
        f.close()
        
        # if the file is completely empty, let's just wipe it!
        if cnt == 0:
            os.remove(path)
           
        if len(list_remove) == 0:
            return []
        else:
            return pd.DataFrame(np.array(list_remove), columns=["Chain", "Resid"]).drop_duplicates()
    
 
    def calculate_pkaANI(self, path):

        code_for_df = os.path.basename(path).split(".")[0]
        pdb_path = path.split(".")[0]

        try:   
            _ = subprocess.run(['pkaani', '-i', path])
        except Exception as e:
            raise Exception("Failed to obtain pkaANI data. %s."%e)

        try:
            log_file = pdb_path + '_pka.log'
            propres = open(log_file)
        except Exception:
            raise Exception('Failed to find pkaANI log file')

        lys_number = list()
        pkas = list()
        chains = list()
        try:
            for line in propres:
                if re.search('^LYS', line):
                    line = line[4:]
                    info = line.split()

                    lys_number.append(int(info[0]))
                    chains.append(info[1])
                    pkas.append(float(info[2]))

            propres.close()
        except Exception as e:
            propres.close()
            raise Exception('Failure parsing %s_pka.log. %s'%(code_for_df, e))
            
        try:
            df = pd.DataFrame({'Resid':lys_number,
                         'Chain': chains,
                         'pkaANI':pkas})

            df.sort_values(by=['pkaANI'], inplace=True)
            df = df.dropna()
            df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e:
            raise Exception('Failed to construct pkaANI dataframe. %s'%e)

        return df
   
 
    def calculate_sasa(self, path):
        '''
        Form small structures which include just the atoms surrounding the lysine of interest,
        and compute the SASA from that.
        '''
        
        try:
            list_of_sasa = list()
            list_of_resid = list()
            list_of_chains = list()
    
            #read PDB file
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            df = M.data
    
            #Find the coordinates and index of all lysine residues in the protein.
            lys_coords, lys_idx = M.atomselect('*', ['LYS'], 'NZ', use_resname=True, get_index=True)
            df = M.data
    
            #Use get_subset...
            #Find the chain and resid of each lysine.
            for entry in lys_idx:
                chain = df.at[entry, 'chain']
                list_of_chains.append(chain)
    
            for entry in lys_idx:
                resid = df.at[entry, 'resid']
                list_of_resid.append(int(resid))
                
            #Find the coordinates and index of every atom in the molecule.
            all_coords, idx = M.atomselect('*','*','*', get_index=True)
    
        except Exception as e:
            raise Exception("%s"%e)
        
        #For each lysine it works out the distance between the lys NZ,
        #and the each atom in the protein.
        for j in range(len(lys_coords)):
            list_close_points = list()
            
            for i in range(len(all_coords)):
                try:
                    x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                    y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                    z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                    distance = np.sqrt(x_dist + y_dist + z_dist)
                    if distance < 15:
                        list_close_points.append(idx[i])
                except:
                    continue
    
            #if the atoms are close to the lys NZ they are included in a small .pdb structure.
            try:
                M.write_pdb('temp_struc.pdb', index=list_close_points, split_struc=False)
    
                S = bb.Molecule()
                S.import_pdb('temp_struc.pdb', include_hetatm=True)
                chain = list_of_chains[j]
                resid = list_of_resid[j]
    
                #SASA is calculated for that lysine in the small molecule.
                pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ"],
                                             use_resname=False, get_index=True)
                x = bb.sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0)
                list_of_sasa.append(x[0])
    
            except:
                print('Error obtaining SASA at index value ' + str(j))
                list_of_sasa.append(None)
                continue
    
        #append results to a df which is given as output
        try:
            df = pd.DataFrame({'Chain': list_of_chains,
                               'Resid': list_of_resid,
                               'sasa': list_of_sasa})
    
        except Exception as e:
            raise Exception('Error obtaining SASA data. %s'%e)
    
        try:
            os.remove('temp_struc.pdb')
            
        except Exception as e:
            print("Error %s"%e)
            print('Failed to remove temporary pdb structure.')
            pass
        
        return df
       
    
    def calculate_depth(self, path):
        
        try:
            M = bb.Molecule(path)
            pos, idx = M.atomselect("*", "*", "NZ", get_index=True)
        except:
            raise Exception("could not find NZ atoms in atomic structure")
        
        try:
            parser = PDBParser()
            structure = parser.get_structure('structure', path)
            surface = get_surface(structure[0])
        except Exception as e:
            raise Exception(f"could not get biopython structure. {e}")
        
        
        results = []
        for i in range(len(pos)):
            chain = M.data.loc[idx[i], ["chain"]].values[0]
            resid = M.data.loc[idx[i], ["resid"]].values[0]
            
            mychain = structure[0][chain]
            myres = mychain[int(resid)]
            
            try:
                #dist = min_dist(pos[i], surface)
                rd = residue_depth(myres, surface)
            except:
                raise Exception(">> failed getting min_dist")
            
            results.append([chain, resid, rd])
         
    
        df_Depth = pd.DataFrame(results, columns=["Chain", "Resid", "depth"])
    
        return df_Depth
    
    

if __name__ == "__main__":

    f1 = "Demo{os.sep}curated{os.sep}1M2E-alt-1.pdb"

    from uniprot import Uniprot
    from protein import PDB

    print("Scanning UNIPROT...")
    UP = Uniprot()
    UP.get_protein_data("P0CG47")

    UP.df = UP.df.iloc[3:10]

    if True:
        print("Gathering proteins")
        pdb = PDB(gap=10)
        pdb.gather_proteins(UP.df)
        pdb.save_state()

    else:
        pdb = PDB(gap=10)
        pdb.load_state('Demo/proteins.csv')

    print(pdb.df)

    print("Measuring...")
    M = Measure(pdb.df)
    M.measure_dataframe()






