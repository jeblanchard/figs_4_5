#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""A module that makes downloading and using datasets from the Gene Expression Omnibus (GEO) easy.

Author: Joshua Blanchard, jeblanchard@berkeley.edu

Currently has limited functionality."""


# In[2]:


import pandas as pd
import os
import numpy as np
import re as re
import xmltodict
import ftplib
import tarfile
import pickle
import gzip


# In[3]:


__base_path = "data"


# In[4]:


__hidden_path = ".aidp_files"


# In[5]:


if not os.path.exists(__base_path):
    os.mkdir(__base_path)
    
if not os.path.exists(__hidden_path):
    os.mkdir(__hidden_path)


# In[6]:


def __family_path(GSE_family):
    
    """Exists to build paths to the family's directory.
    
    I downloaded the files using the download() function in conjuction with the extract() function."""
        
    return os.path.join(__base_path, GSE_family + "/")


# In[7]:


def __clean(file_path):
    
    """Cleans a given .txt file.
    
    Returns a dictionary:
    
    "site": first column
    "measurement": second column
    "bad_rows": list of all the invalid rows"""

    valid_rows = []
    not_valid_rows = []
    file = open(file_path, 'r')
    
    for line in file:
        
#         checks for only the first two columns
        line_match = re.match(r"\S+\t\S+", line)
        if line_match:
            valid_rows.append(line_match.group(0))
        else:
            not_valid_rows.append(line)
        
    file.close()
    
#     now let's split our valid_rows list into two lists, one for each column
    col_1 = []
    col_2 = []
    for row in valid_rows:
        row_match = re.match(r"(\S+)\t(\S+)", row)
        col_1.append(row_match.group(1))
        col_2.append(row_match.group(2))
        
    return {"col_1":col_1, "col_2":col_2, "bad_rows":not_valid_rows}


# In[8]:


def __load_file(file_directory):
#     clean data file
#     convert to a dataframe object and return

    """Given a file name will output a corresponding pandas.DataFrame object."""
    
    try:
        clean_dict = __clean(file_directory)
    except PermissionError:
        print("You likely inputted the path of a directory, not a file.")
    
#     return pd.DataFrame({"site": clean_dict["col_1"], "measurement": clean_dict["col_2"]})
    return pd.DataFrame(data= clean_dict["col_2"], index= clean_dict["col_1"], columns= ["measurement"])


# In[9]:


def family_dict(GSE_family):
    
    """*** DEPRECATED *** 
    
    Given a family ID, will output a dictionary. Keys will be the sample IDs, values will be the corresponding
    pandas.DataFrame object."""

#     let's check if we already have this dictionary saved
    if not os.path.exists("./" + __hidden_path):
        os.mkdir(__hidden_path)
    
    dict_path = __hidden_path + "/" + GSE_family + "_dict"    
    if not os.path.exists(dict_path):
    
        family_directory = __family_path(GSE_family)
        total_list = os.listdir(family_directory)
        valid_files = []

        for file_name in total_list:
            match = re.match(r"GSM", file_name)
            if match:
                valid_files.append(file_name)

        family_dict = {}
        for file_name in valid_files:
            file_df = __load_file(os.path.join(family_directory, file_name))
            sample_id = re.match(r"GSM\d+", file_name).group(0)
            family_dict[sample_id] = file_df
            
        dict_file = open(dict_path, 'wb')
        pickle.dump(family_dict, dict_file)
        dict_file.close()
        
    else:
        dict_file = open(dict_path, 'rb')
        family_dict = pickle.load(dict_file)
        dict_file.close()
        
        
    return family_dict


# In[10]:


def __matrix_helper(file_path, use= "series()"):
    
    """Returns a tuple containing the start line for reading (0) and the number of rows to read (1).
    
    Possible values for use: "series()", "info()"
    """
    
    file = open(file_path, mode= 'r', errors= 'replace')
    
    line_num = 0
    if use == "series()":
        
        while True:
            line = file.readline()
            
            if "!series_matrix_table_begin" in line:
                start_row = line_num
            elif "!series_matrix_table_end" in line:
                end_row = line_num - 1
            elif line == "":
                break

            line_num += 1
        
    else:
                
        while True:
            line = file.readline()

            if "!Sample_title" in line:
                start_row = line_num
            elif "!series_matrix_table_begin" in line:
                end_row = line_num - 1
            elif line == "":
                break

            line_num += 1

    num_rows = end_row - start_row - 1
        
    return start_row, num_rows


# In[11]:


def __matrix_to_df(file_path, use= "series()", GSE= ""):
    
    """Returns the pandas.dataframe corresponding to the series_matrix file.
    
    Possible values for use: "series()", "info()"
    """

    if use == "series()":
        
        series_path = os.path.join(__hidden_path + "/" + GSE + "/" + "series_df")
    
        if os.path.exists(series_path):

            df_file = open(series_path, "rb")
            df = pickle.load(df_file)
            df_file.close()
        
        else:
            
            start_row, num_rows = __matrix_helper(file_path, use)
            df = pd.read_csv(file_path, header= start_row, sep= "\t", low_memory= False, nrows= num_rows)
            df.set_index("ID_REF", inplace= True)
            
            GSE_dir = os.path.join(__hidden_path + "/" + GSE)
            if not os.path.exists(GSE_dir):
                os.mkdir(GSE_dir)
            
            df_file = open(series_path, 'wb')
            pickle.dump(df, df_file)
            df_file.close()
            
    elif use == "info()":
        
        info_path = os.path.join(__hidden_path + "/" + GSE + "/" + "info_df")
    
        if os.path.exists(info_path):

            info_df_file = open(info_path, "rb")
            df = pickle.load(info_df_file)
            info_df_file.close()
    
        else:
            
            start_row, num_rows = __matrix_helper(file_path, use)
            df = pd.read_csv(file_path, header= start_row, sep= "\t", low_memory= False, nrows= num_rows)
            
            df.set_index("!Sample_geo_accession", inplace= True)
            df = df.loc["!Sample_characteristics_ch1"]
            
            GSE_dir = os.path.join(__hidden_path + "/" + GSE)
            if not os.path.exists(GSE_dir):
                os.mkdir(GSE_dir)
            
            df_file = open(info_path, 'wb')
            pickle.dump(df, df_file)
            df_file.close()

    return df


# In[12]:


def info(GSE_family, sample_id, info= "age"):
    
    """Given a GSE family and the ID of the wanted sample will return the desired information of the sample.
    
    Possible values for the info parameter:
    
    "age": the unit of the outputted age will always be years.
    "brca1": the BRCA1 mutation status for GSE57285
    "arthritis": the arthritis status for GSE42861
    "crohns": the Crohn's disease status for GSE32148"""
#     create dataframe of sample characteristics portion of the series_matrix file
# read desired info from this dataframe. will have to use regex to find the right information
# set up some sort of persistance of this dataframe for faster retrieval of information

#             match = re.search(r"(a|A)(g|G)(e|E)", tag)

    info_df = __matrix_to_df("./" + __base_path + "/" + GSE_family + "/" + GSE_family + "_series_matrix.txt", use= "info()", GSE= GSE_family)    
    info_series = info_df.loc[:, sample_id]
        
    
    if info == "age":
        
        for row in info_series:
            match_age = re.search(r"(a|A)(g|G)(e|E)", row)
            if match_age:
        #         let's grab the age
                match_age_dig = re.search(r"\d+\.?\d*", row)
                if match_age_dig:
                    age_float = float(match_age_dig.group())
    #             to account for "newborn" instead of an actual number
                else:
                    age_float = 0.0

        #         now let's find the unit of age
                match_years = re.search(r"(y|Y)(e|E)(a|A)(r|R)(s|S)*", row)
                match_months = re.search(r"(m|M)(o|O)(n|N)(t|T)(h|H)(s|S)*", row)
                match_days = re.search(r"(d|D)(a|A)(y|Y)(s|S)*", row)
                match_hours = re.search(r"(h|H)(o|O)(u|U)(r|R)(s|S)*", row)

        #         we always return the age in units of years
                if match_years:
                    return age_float
                elif match_months:
                    return age_float / 12
                elif match_days:
                    return age_float / 365
                elif match_hours:
                    return age_float / 8760
                else:
        #             if no unit is given we'll default to years
                    return age_float
        
    elif info == "brca1":
            
        for row in info_series:
            match_status = re.search(r"brca1 mutation status", row)
            if match_status:
        #         let's grab the brca1 status
                match_status_dig = re.search(r"(brca1 mutation status: )(\d)", row)
                if match_status_dig:
                    status_int = int(match_status_dig.group(2))

                    return status_int
    
    elif info == "arthritis":
        
        for row in info_series:
            match_arth = re.search(r"disease state", row)
            if match_arth:
        #         let's grab the arthritis status
                match_arth_status = re.search(r"(disease state: )(\S+)", row)
                if match_arth_status.group(2) == "rheumatoid":
                    return "rheumatoid arthritis"
                else:
                    return "normal"
                
    else:
        
        for row in info_series:
            match_crohns = re.search(r"disease state", row)
            if match_crohns:
        #         let's grab the Crohn's status
                match_crohns_status = re.search(r"(disease state: )(\S+)", row)
                return match_crohns_status.group(2)
        


# In[13]:


def __ID_to_int(GSE_ID):
    
    """Given some GSE ID will return the corresponding integer as an int.
    
    Example:
    
    __ID_to_int("GSE41037") will return 41037."""
    
    match = re.search(r"\d+", GSE_ID)
    
    if not match:
        print(GSE_ID)
    
    return int(match.group(0))


# In[14]:


def __sub_directory(GSE_ID):
    
    """Given a GSE ID will return the corresponding sub-directory."""
    
    gse_int = __ID_to_int(GSE_ID)
    
    if gse_int <= 171:
        ret_str = "GSE" + str(gse_int) + "nnn"
    else:
        first_3_dig = int(str(gse_int)[0:3])
        if first_3_dig <= 171:
            ret_str = "GSE" + str(first_3_dig) + "nnn"
        else:
            first_2_dig_str = str(gse_int)[0:2]
            ret_str = "GSE" + first_2_dig_str + "nnn"
            
    return ret_str


# In[15]:


def download(GSE_family_list, file_type= "miniml"):
    
    """Will download family .tgz files to the following directory: "./data/"
    
    Possible values for file_type: "miniml", "series_matrix"
    """
    
    if (type(GSE_family_list) == str):
        GSE_family_list = [GSE_family_list]

    url = "ftp.ncbi.nlm.nih.gov"
    ftp = ftplib.FTP(url)
    ftp.login()
    
    if not os.path.exists(__base_path):
        os.mkdir(__base_path)
    
    for GSE_family in GSE_family_list:
        
        if file_type == "miniml":
        
            ftp.cwd("/geo/series/" + __sub_directory(GSE_family) + "/" + GSE_family + "/miniml/")
            filename = GSE_family + "_family.xml.tgz"
            
            if not (os.path.exists(__base_path + "/" + filename) or os.path.exists(__base_path + "/" + GSE_family)):
                local_file = open(__base_path + "/" + filename, 'wb')
                ftp.retrbinary('RETR ' + filename, local_file.write, blocksize= 16_384)
                local_file.close()
            
        else:
            
            ftp.cwd("/geo/series/" + __sub_directory(GSE_family) + "/" + GSE_family + "/matrix/")
            filename = GSE_family + "_series_matrix.txt.gz"
        
            if not (os.path.exists(__base_path + "/" + filename) or os.path.exists(__base_path + "/" + GSE_family + "/" + GSE_family + "_series_matrix.txt")):
                local_file = open(__base_path + "/" + filename, 'wb')
                ftp.retrbinary('RETR ' + filename, local_file.write, blocksize= 16_384)
                local_file.close()

    ftp.quit()
    


# In[16]:


def __check_beta(series_df):
    
    """Returns True if Beta values are recorded, returns False if M values are recorded."""


# In[17]:


def series(GSE):
    
    """Returns a pandas.DataFrame representative of the series' data."""
    
#     download(GSE, file_type= "series_matrix")
#     __extract()
    
    file_path = os.path.join(__base_path, GSE, GSE + "_series_matrix.txt")
    
    return __matrix_to_df(file_path, GSE= GSE)


# In[18]:


def __xml_path(GSE_family):
    
    """Exists to build a path to the family's .xml.
    
    I downloaded these families using the download() function in conjuction with the extract() function."""
    
    return os.path.join(__base_path, GSE_family + "/" + GSE_family + "_family.xml")
    


# In[19]:


def __xml_to_dict(xml_path):
    
    """Exists to convert a .xml file at xml_path to a dictionary."""
    
    family_file = open(xml_path,'r+b')
    family_dict = xmltodict.parse(family_file)
    family_file.close()
    
    return family_dict


# In[20]:


def __dict_index(GSE_family, sample_id):
    
    """Gives the index of where in the dictionary the sample's information is."""
    
    return __sample_indices(GSE_family)[sample_id]


# In[21]:


def __sample_indices(GSE_family):
    
    """Exists to return the indices of each sample's information within the associated family dictionary. Returns a dictionary
    with keys equal to the sample ID ("GSM***") and values equal to the index of that sample's information within the family
    dictionary."""
    
    family_dir = __family_path(GSE_family)
    file_list = os.listdir(family_dir)
    
    filtered_list = []
    for file_name in file_list:
        sample_match = re.match(r"GSM", file_name)
        if sample_match:
            filtered_list.append(file_name)
            
    for i in np.arange(len(filtered_list)):
        filtered_list[i] = re.match(r"GSM\d+", filtered_list[i]).group(0)
        
    index_dict = {}
    index = 0

    for sample_id in filtered_list:
        index_dict[sample_id] = index
        index += 1
        
    return index_dict


# In[22]:


def __extract():
    
    """Will extract files from all downloaded family .tgz files to a respective directory: ./data/GSE***/
    
    This will also delete the .tgz and .gz files."""

    file_list = os.listdir(__base_path)
    tgz_list = []
    gz_list = []
    
    for filename in file_list:
        match_tgz = re.search(r"\.tgz", filename)
        match_gz = re.search(r"\.gz", filename)
        
        if match_tgz:
            tgz_list.append(filename)
        elif match_gz:
            gz_list.append(filename)
                
    flag = False
    
    for filename in tgz_list:
        full_path = __base_path + "/" + filename
        family_id = re.search(r"GSE\d+", filename).group(0)
                
        file = tarfile.open(full_path)
        
        try:
            out_path = "./" + __base_path + "/" + family_id
            file.extractall(out_path)
        except:
#             let's end the content extraction of the file. will fix later.
            file.close()
    
            flag = True
            print("Something weird happened while extracting from the " + family_id + " compressed file. Ended extraction early for " + family_id + ".")
            
        if not flag:
            file.close()
            os.remove(full_path)
            
    
    for filename in gz_list:
        full_path = "./" + __base_path + "/" + filename
        family_id = re.search(r"GSE\d+", filename).group(0)
                
        compressed_file = open(full_path, 'rb')
        compressed_file_contents = compressed_file.read()
        compressed_file.close()
                
        split_compressed_bytes = compressed_file_contents.split()
        contents_str = ''
        
        
        for compressed_bytes in split_compressed_bytes:
            contents_bytes = gzip.decompress(compressed_bytes)
            contents_str.append(contents_bytes.decode(errors= "replace"))
        
        family_dir = os.path.join(__base_path, family_id)
        if not os.path.exists(family_dir):
            os.mkdir(family_dir)
        
        out_path = "./" + __base_path + "/" + family_id + "/" + family_id + "_series_matrix.txt"

        file = open(out_path, mode= 'w', encoding= "utf-8")
        
        file.write(contents_str)
        file.close()
    
        if not flag:
            file.close()
            os.remove(full_path)
    

