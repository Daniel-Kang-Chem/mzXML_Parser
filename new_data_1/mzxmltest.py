from pyteomics import mzxml, pepxml
import numpy as np
import os
import time
import glob as glob
import matplotlib.pyplot as plt
MAX_SCANS = 400  # maximum number of mzxml scans to look at once (must be greater than 50)
RETENTION_TIME_DELTA = 30  # retention time window to look in for summation integration
DELTA_MZ = 0.005  # adjusts the resolution of the mz array (filling in data points with interpolated values)
PEPTIDES = ["sp|P00359|G3P3_YEAST"]

def find_mzxml_pepxml():
    """
    function to find pepxml and mzxml files within the folder
    :return: tuple of two elements, first one which is the first mzxml file found and the second which is the first pep
    xml file found
    """
    mzxml_file_name = glob.glob("*.mz*")
    pepxml_file_name = glob.glob("*.pep.xml")
    return (mzxml_file_name,pepxml_file_name)


def list_summation_integration(scan_list, mass, charge, minrt, maxrt, peptide_name, dmz, label=None, min_mz=None, max_mz=None):
    """
    prints out a scatter plot of summation integration as per Davis's requirements
    :param scan_list: list of scans
    :param mass: mass found by the pepxml file
    :param charge: assumed charge from pepxml file
    :param minrt: minimum retention time to plot from
    :param maxrt: maximum retention time to plot to
    :param peptide_name: name of the peptide being graphed
    :param dmz: delta mz used to create y axis (default=.05)
    :param min_mz: minimum mz to graph from
    :param max_mz: maximum mz to graph to
    :param label: "SILAC-labeling"
    :return: prints out a scatter plot using parameters and Davis's requirements
    """

    if label == "silac-labeling":
#TODO: look at sequence, +8 lysine and +10 arginine
        silac_lower = int((mass+charge)/charge-5/charge) # always 5/charge for lower shift
        silac_upper = int((mass+charge)/charge+25/charge)
        return list_summation_integration(scan_list, mass, charge, minrt, maxrt, peptide_name, dmz, None, silac_lower, silac_upper)
    elif label == "15n-labeling":
#TODO: make a dictionary of amino acids and how many nitrogens it has, use to calculate m/z shift by looking at peptide sequence
        return list_summation_integration()
    elif label == "13c-labeling":
#TODO: same as above, but for carbon
        return list_summation_integration()

    mzarr = np.arange(float(min_mz),float(max_mz),dmz)
    mzlen = len(mzarr)
    intensity = np.zeros(mzlen)
    start_index = 0
    end_index = len(scan_list) - 1
    search_index = int(end_index/2)

    while scan_list[search_index]['retentionTime'] * 60 > minrt:
        end_index = search_index
        search_index = int((start_index+end_index)/2)
    while scan_list[search_index]['retentionTime'] * 60 < minrt:
        start_index = search_index
        search_index = int((start_index+end_index)/2)
    while scan_list[search_index]['retentionTime'] * 60 > minrt:
        search_index -= 1
    while scan_list[search_index]['retentionTime'] * 60 < minrt:
        search_index += 1

    i = int(search_index)

    while scan_list[i]['retentionTime'] * 60 <= maxrt:

        #interpIntensArr = np.interp(mzarr, scan_list[i]['m/z array'], scan_list[i]['intensity array'])
        interp_intens_arr = interpolate_scan(scan_list[i], mass, charge, dmz, "silac-labeling")
        for x in range(0,mzlen-1):
            intensity[x] += interp_intens_arr[x]
        i += 1

    filename = 'scans/test_data_' + str(peptide_name) + '_' + str(mass) + '.txt'

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename,'w+') as f:
        for j in range(len(intensity)):
            f.write(str(mzarr[j]) + " " + str(intensity[j]) + "\n")

    print("finished with:",str(peptide_name))


def interpolate_scan(scan, mass, charge, dmz, label=None, min_mz=None, max_mz=None):
    """
    returns a scan dictionary with its mz/intensity array interpolated to fit DELTA_MZ
    :param scan: dictionary object given by mzxml containing mz and intensity array
    :return: scan dictionary object with mz and intensity array interpolated to fit DELTA_MZ increments of mz
    """
    if label == "silac-labeling":
        silac_lower = int((mass+charge)/charge-5/charge)
        silac_upper = int((mass+charge)/charge+25/charge)
        return interpolate_scan(scan, mass, charge, dmz, None, silac_lower, silac_upper)
    elif label == "15n-labeling":
        return 1
    elif label == "13c-labeling":
        return 1
    mzarr = np.arange(float(min_mz), float(max_mz), dmz)
    interpIntensArr = np.interp(mzarr, scan['m/z array'], scan['intensity array'])
    scan['m/z array'] = mzarr
    scan['intensity array'] = interpIntensArr
    return scan

def find_first_hit(iter,minrt):
    """
    returns an iterator that starts with the first scan after the minimum retention_time
    :param iter: mzxmlfile iterator
    :param minrt: minimum retention_time (from first pepxml hit)
    :return: iterator at first scan after the minimum retention time
    """
    scan = iter.next()
    while scan['retentionTime']*60 < minrt:
        scan = iter.next()
    return iter


def run_program(mzxml_iterator,pepxml_iterator,scan_list=[]):
    try:
        pep_scan = pepxml_iterator.next()
    except:
        print("finished with pepxml")
        return time.time()

    mzxml_iterator = find_first_hit(mzxml_iterator, pep_scan['retention_time_sec'] - RETENTION_TIME_DELTA)  # skip scans that will not be used.

    scan_list_len = len(scan_list)
    for i in range(MAX_SCANS-scan_list_len):
        try:
            scan = mzxml_iterator.next()
        except:
            print("finished mzxml file")
            return time.time()
        if i == MAX_SCANS - 50 - scan_list_len:
            ret_list = [scan]
        elif i > MAX_SCANS - 50 - scan_list_len:
            ret_list.append(scan)
        scan_list.append(scan)

    while pep_scan['retention_time_sec'] < scan_list[-1]['retentionTime'] * 60 - 30:
        peptide_name = pep_scan['search_hit'][0]['proteins'][0]['protein']
        if "RAND" not in peptide_name and peptide_name in PEPTIDES:  # look for peptides

            print("Peptide: " + pep_scan['search_hit'][0]['peptide'])

            id_time = pep_scan['retention_time_sec']
            precursor_mass = pep_scan['precursor_neutral_mass']
            assumed_charge = pep_scan['assumed_charge']
            # TODO: Find smart retention time to integrate across (replacing rt delta parameter)
            # TODO: We have a scan list, the pepxml rt time, the pepxml assumed charge, pepxml neutral mass
            mz = (precursor_mass + assumed_charge) / assumed_charge
            mz = round(mz, 3) 
            t_mz = mz * 1000
            t_mz += 5 - t_mz % 5

            mz = t_mz / 1000

            int_rt_arr = []
            rt_arr = []

            # interpolate the scan_list
            for index, value in enumerate(scan_list):

                scan_list[index] = interpolate_scan(scan_list[index], precursor_mass, assumed_charge, DELTA_MZ, "silac-labeling")

                print(mz)
                # find the index to pull from intensity array by finding the mz in m/z array
                for index2, scan_mz in enumerate(scan_list[index]['m/z array']):
                    print(scan_mz)
                    if round(scan_mz,3) == round(mz,3):
                        index_found = index2
                     
                int_rt_arr.append(scan_list[index]['intensity array'][index_found])
                rt_arr.append(scan_list[index]['retentionTime'])

            plt.plot(rt_arr, int_rt_arr)
            plt.show()

            start_rt = pep_scan['retention_time_sec']

            mz_found = int(
                (pep_scan['precursor_neutral_mass'] + pep_scan['assumed_charge']) / pep_scan['assumed_charge'])

            list_summation_integration(scan_list, pep_scan['precursor_neutral_mass'], pep_scan['assumed_charge'],
                                       start_rt - RETENTION_TIME_DELTA, start_rt + RETENTION_TIME_DELTA,
                                       peptide_name,
                                       DELTA_MZ, "silac-labeling")  # HARD CODED SILAC LABELING



            data_file_name = 'scans/test_data_' + str(peptide_name) + '_' + str(pep_scan['precursor_neutral_mass']) + '.txt'

            with open('misc/NPULSE.batch', 'a') as f:
                f.write(str(pep_scan['search_hit'][0]['peptide']) + " " + str(pep_scan['assumed_charge']) + " " + str(data_file_name) + "\n")

        try:
            pep_scan = pepxml_iterator.next()
        except:
            print("finished with pepxml file")
            return time.time()
    del scan_list[:]
    del scan_list
    return run_program(mzxml_iterator,pepxml_iterator,ret_list)

results = find_mzxml_pepxml()
print("Found mzxml file:", results[0][0], "and pepxml file", results[1][0])

mzxmlfilename = results[0][0]
print("accessing mzxml file:", mzxmlfilename)
mzxml_it = mzxml.MzXML(mzxmlfilename)

pepxmlfilename = results[1][0]
print('accessing pepxml file:', pepxmlfilename)
try:
    pepxml_it = pepxml.PepXML(pepxmlfilename)
except:
    pepxml_it = pepxml.PepXML(pepxmlfilename)
os.makedirs(os.path.dirname('misc/NPULSE.batch'), exist_ok=True)  # create NPULSE file
t1 = time.time()
t2 = run_program(mzxml_it, pepxml_it)
print(t1,t2)
print(t2 - t1)
dt = t2 - t1
print("time taken:",str(dt))
