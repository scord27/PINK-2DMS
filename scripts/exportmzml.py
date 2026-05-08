# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 12:01:38 2024

@author: lwg513
"""

import sys
from tkinter import filedialog
sys.path.append("pinkpeak")
from psims.mzml.writer import MzMLWriter
import numpy as np

class exportasmzml(): #Klein, J. A., & Zaia, J. (2018). psims - A declarative writer for mzML and mzIdentML for Python. Molecular & Cellular Proteomics, mcp.RP118.001070. https://doi.org/10.1074/mcp.RP118.001070
    def write_to():
        mzMLsaveas=filedialog.asksaveasfile(mode='wb',defaultextension=".mzML",filetypes=[("mzML Files", "*.mzML")])
        return mzMLsaveas
    def writemzml(mzMLsaveas,data,source_type,analyser_type,Polarity,fragmentation_technique,fragmentation_power,Fragmentation_duration,MS_level,ExD_bias, ExD_lens,precursor_int):
        writer=MzMLWriter(mzMLsaveas)
        with writer:
            writer.controlled_vocabularies()
            writer.file_description(["MSn spectrum","centroid spectrum"])
            writer.software_list([ {"id": "psims-writer","version": "0.1.2","params": ["python-psims",]}])
            source = writer.Source(1, [source_type, "electrospray inlet"])
            analyzer = writer.Analyzer(2, [analyser_type])
            detector = writer.Detector(3, ["inductive detector"])
            config = writer.InstrumentConfiguration(id="IC1", component_list=[source, analyzer, detector],params=["Bruker Daltonics instrument model"])
            writer.instrument_configuration_list([config])
            methods = []
            methods.append(writer.ProcessingMethod(order=1, software_reference="psims-writer", params=tuple(["peak picking", "Conversion to mzML"])))
            processing = writer.DataProcessing(methods, id='DP1')
            writer.data_processing_list([processing])
            with writer.run(id=1, instrument_configuration='IC1'):
                with writer.spectrum_list(count=len(data)):
                    i = 0
                    while i < len(data):
                        msn_mzs=[]
                        msn_intensities=[]
                        for j in range (1,len(data[i])):
                            if len(data[i])<2:
                                continue
                            else:
                                msn_mzs.append(data[i][j][0])
                                msn_intensities.append(data[i][j][2])
                        msn_intensities=np.array(msn_intensities)
                        msn_spectrum_id = "index=%d" % i
                        scan_time = 0.23 * i
                        precursor = writer.precursor_builder()
                        si = precursor.selected_ion(mz=str(precursor_int[0][i]),intensity=str(precursor_int[1][i]))
                        si.charge = 1
                        act = precursor.activation()
                        if fragmentation_technique=="IRMPD":
                            act.add_param(fragmentation_technique)
                            act.add_param({"percentage": fragmentation_power})
                            act.add_param({"wavelength": 10.1, "unitName":"micrometer"})
                            act.add_param({"laser power": 25, "unitName":"watt"})
                            act.add_param({"pulse duration": Fragmentation_duration, "unitName":"second"})
                        elif fragmentation_technique=="ExD":
                            act.add_param(fragmentation_technique)
                            act.add_param({"Bias": 10, "unitName":"volt"})
                            act.add_param({"Lens": 10, "unitName":"volt"})
                            act.add_param({"pulse duration": Fragmentation_duration, "unitName":"second"})
                        elif fragmentation_technique=="UVPD":
                            act.add_param(fragmentation_technique)
                            act.add_param({"percentage": fragmentation_power})
                            act.add_param({"wavelength": 10.1, "unitName":"micrometer"})
                            act.add_param({"laser power": 25, "unitName":"watt"})
                            act.add_param({"pulse duration": Fragmentation_duration, "unitName":"second"})
                        else:
                            print("please specify a fragmentation technique")
                            return
                        writer.write_spectrum(msn_mzs, msn_intensities, id=msn_spectrum_id,polarity=Polarity, centroided=True,precursor_information=precursor,params=[{"ms level": MS_level}, {"total ion current": msn_intensities.sum()}])
                        i=i+1
        return "finished exporting"
    def writesinglemzml(mzMLsaveas,data,source_type,analyser_type,Polarity,fragmentation_technique,fragmentation_power,Fragmentation_duration,MS_level,ExD_bias, ExD_lens):
        writer=MzMLWriter(mzMLsaveas)
        with writer:
            writer.controlled_vocabularies()
            writer.file_description(["MSn spectrum","centroid spectrum"])
            writer.software_list([ {"id": "psims-writer","version": "0.1.2","params": ["python-psims",]}])
            source = writer.Source(1, [source_type, "electrospray inlet"])
            analyzer = writer.Analyzer(2, [analyser_type])
            detector = writer.Detector(3, ["inductive detector"])
            config = writer.InstrumentConfiguration(id="IC1", component_list=[source, analyzer, detector],params=["Bruker Daltonics instrument model"])
            writer.instrument_configuration_list([config])
            methods = []
            methods.append(writer.ProcessingMethod(order=1, software_reference="psims-writer", params=tuple(["peak picking", "Conversion to mzML"])))
            processing = writer.DataProcessing(methods, id='DP1')
            writer.data_processing_list([processing])
            with writer.run(id=1, instrument_configuration='IC1'):
                with writer.spectrum_list(count=1):
                    i = 0
                    while i < 1:
                        msn_mzs=[]
                        msn_intensities=[]
                        for j in range (1,len(data[0])):
                            msn_mzs.append(data[0][j])
                            msn_intensities.append(data[2][j])
                        msn_intensities=np.array(msn_intensities)
                        msn_spectrum_id = "index=%d" % i
                        scan_time = 0.23 * i
                        precursor = writer.precursor_builder()
                        si = precursor.selected_ion(mz=str(data[0][i]).strip("[]"))
                        si.charge = 1
                        act = precursor.activation()
                        if fragmentation_technique=="IRMPD":
                            act.add_param(fragmentation_technique)
                            act.add_param({"percentage": fragmentation_power})
                            act.add_param({"wavelength": 10.1, "unitName":"micrometer"})
                            act.add_param({"laser power": 25, "unitName":"watt"})
                            act.add_param({"pulse duration": Fragmentation_duration, "unitName":"second"})
                        elif fragmentation_technique=="ExD":
                            act.add_param(fragmentation_technique)
                            act.add_param({"Bias": 10, "unitName":"volt"})
                            act.add_param({"Lens": 10, "unitName":"volt"})
                            act.add_param({"pulse duration": Fragmentation_duration, "unitName":"second"})
                        elif fragmentation_technique=="UVPD":
                            act.add_param(fragmentation_technique)
                            act.add_param({"percentage": fragmentation_power})
                            act.add_param({"wavelength": 10.1, "unitName":"micrometer"})
                            act.add_param({"laser power": 25, "unitName":"watt"})
                            act.add_param({"pulse duration": Fragmentation_duration, "unitName":"second"})
                        else:
                            print("please specify a fragmentation technique")
                            return
                        writer.write_spectrum(msn_mzs, msn_intensities, id=msn_spectrum_id,polarity=Polarity, centroided=True,precursor_information=precursor,params=[{"ms level": MS_level}, {"total ion current": msn_intensities.sum()}])
                        i=i+1
        return "finished exporting"
        