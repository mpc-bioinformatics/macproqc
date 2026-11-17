import streamlit as st
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import os
import subprocess






def select_folder():
   root = tk.Tk()
   root.withdraw()
   folder_path = filedialog.askdirectory(master=root)
   root.destroy()
   return folder_path

def select_file():
   root = tk.Tk()
   root.withdraw()
   folder_path = filedialog.askopenfilename(master=root)
   root.destroy()
   return folder_path

st.title("Next-QC-Flow - Simple GUI")
st.subheader("The next quality control workflow, recreated in Nextflow.")
st.divider()

### MAIN PARAMETERS
# Select Input Folders and Files
st.subheader("Input-Parameters (REQUIRED)")
selected_RAW_folder = st.session_state.get("selected_RAW_folder", None)
folder_select_RAWS = st.button("Select *.RAW/*.d-files")
if folder_select_RAWS:
  selected_RAW_folder = select_folder()
  st.session_state.selected_RAW_folder = selected_RAW_folder
if selected_RAW_folder:
   st.write("Selected folder path:", selected_RAW_folder)

selected_FASTA_file = st.session_state.get("selected_FASTA_file", None)
file_select_FASTA = st.button("Select FASTA-File")
if file_select_FASTA:
  selected_FASTA_file = select_file()
  st.session_state.selected_FASTA_file = selected_FASTA_file
if selected_FASTA_file:
   st.write("Selected FASTA:", selected_FASTA_file)

selected_COMET_file = st.session_state.get("selected_COMET_file", None)
file_select_COMET = st.button("Select Comet-Search Parameters file")
if file_select_COMET:
  selected_COMET_file = select_file()
  st.session_state.selected_COMET_file = selected_COMET_file
if selected_COMET_file:
   st.write("Selected Search parameters", selected_COMET_file)

selected_OUTPUT_folder = st.session_state.get("selected_OUTPUT_folder", None)
folder_select_OUTPUT = st.button("Select Output-Folder")
if folder_select_OUTPUT:
  selected_OUTPUT_folder = select_folder()
  st.session_state.selected_OUTPUT_folder = selected_OUTPUT_folder
if selected_OUTPUT_folder:
   st.write("Selected Output-folder:", selected_OUTPUT_folder)
st.divider()


### ADDITIONAL PARAMETERS
st.subheader("Additional Parameters")
st.write("Include here additional parameters. E.G.:  Set \"-ccff_header_in_raws\" as a Parameter, and \"-thtp 'Ambient temp. (°C)'\" as a value to extract the ambient temperature from Thermo-RAWS")
st.write("Here, every parameter can be configured for the QC-Workflow!")
df = pd.DataFrame(columns=["Parameter", "Value"])
st.data_editor(df, num_rows="dynamic", key="additional_params")
st.divider()

### Nextflow specific parameters
st.subheader("Nextflow-specific Parameters")
with st.expander("Only Modify (the temp folder) if you know what you are doing"):
    selected_TEMP_folder = st.session_state.get("selected_TEMP_folder", os.path.dirname(__file__) + os.sep + "work")
    folder_select_TEMP = st.button("Select TEMP-Folder")
    if folder_select_TEMP:
        selected_TEMP_folder = select_folder()
        st.session_state.selected_TEMP_folder = selected_TEMP_folder
    if selected_TEMP_folder:
        st.write("Selected TEMP-folder:", selected_TEMP_folder)
    # This directory corresponds to the env var: NXF_TEMP
    with_timeline = st.checkbox("Export Nextflow-Timeline  Report", value=True)  # -with-timeline?
    with_execution = st.checkbox("Export Nextflow Execution Report", value=True)  # -with-report?
st.divider()

### Run Quality Control Workflow!
st.subheader("Execute Next-QC-Flow")
run_qc = st.button("Start the QC in a console!")
if run_qc:
    runnable = True
    # CHECK if everything is set
    if selected_RAW_folder is None or not os.path.exists(selected_RAW_folder):
        st.warning("RAW-File-Folder does not exists!")
        runnable = False
    if selected_FASTA_file is None or not os.path.exists(selected_FASTA_file):
        st.warning("FASTA-file does not exists!")
        runnable = False
    if selected_COMET_file is None or not os.path.exists(selected_COMET_file):
        st.warning("Comet Parameters File does not exists!")
        runnable = False
    if selected_OUTPUT_folder is None or not os.path.exists(selected_OUTPUT_folder):
        _ = os.makedirs(selected_OUTPUT_folder) if selected_OUTPUT_folder is not None else None
    if selected_TEMP_folder is None or not os.path.exists(selected_TEMP_folder):
        _ = os.makedirs(selected_TEMP_folder) if selected_TEMP_folder is not None else None

    if runnable:
        # set initial process string
        qc_exec_str = "pipenv run nextflow run main.nf" + \
            " --main_raw_spectra_folder " + selected_RAW_folder + \
            " --main_fasta_file " + selected_FASTA_file + \
            " --main_comet_params " + selected_COMET_file + \
            " --main_outdir " + selected_OUTPUT_folder + \
            " -work-dir " + selected_TEMP_folder + " "

        # Append additional parameters
        for x in st.session_state["additional_params"]["added_rows"]:
            qc_exec_str += x["Parameter"] + " " + x["Value"] + " "

        if with_timeline:
            timelinefile = selected_OUTPUT_folder + os.sep + "timeline.html"
            qc_exec_str += "-with-timeline " + timelinefile + " "
            _ = os.remove(timelinefile) if os.path.exists(timelinefile) else None
        if with_execution:
            reportfile = selected_OUTPUT_folder + os.sep + "report.html"
            qc_exec_str += "-with-report " + reportfile + " "
            _ = os.remove(reportfile) if os.path.exists(reportfile) else None



        st.write(qc_exec_str)
        st.write("A windows was opened, executing this qc!")
    
        subprocess.run(["/bin/konsole", "--hold", "--separate", "-e", qc_exec_str])