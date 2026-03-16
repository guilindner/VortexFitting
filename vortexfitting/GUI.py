#!/usr/bin/env/ python3
"""
All the basics for the Graphical User Interface, with wxPython
"""

import configparser
import argparse
import os
import sys
from functools import partial

from vortexfitting import fitting
from vortexfitting import schemes
from vortexfitting import detection
from vortexfitting import output
from vortexfitting import classes

import convertToASCII
import convertToNC
import generateNetCDF
import wx.adv
import wx
import wx.html
from GUI_utils import *
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

from importlib.resources import files

logo_file = str(files("vortexfitting") / "logo_VF.ico")

def main():
    app = MyApp()
    app.MainLoop()

class MainFrame(wx.Frame):
    def __init__(self, *args, **kw):
        super(MainFrame, self).__init__(*args, **kw)
        icon = wx.Icon(logo_file, wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

        notebook = wx.Notebook(self)

        # Create tabs
        self.main_tab = wx.Panel(notebook)
        self.utilities_tab = UtilitiesTab(notebook)
        self.about_tab = AboutTab(notebook)

        # Add tabs to notebook
        notebook.AddPage(self.main_tab, "Main")
        notebook.AddPage(self.utilities_tab, "Utilities")
        notebook.AddPage(self.about_tab, "About")

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notebook, 1, wx.EXPAND)  # Permet d'étendre le notebook
        self.SetSizer(sizer)
        self.Maximize(True)
        self.Layout()
        self.analysis_parameters = {
            'input_filename': "./data/example_vel_{:06d}.dat ",
            'output_directory': './results/example_temporal_series',
            'file_type': 'PIV Tecplot (.dat)',
            'meanfile': './data/example_mean.dat',
            'scheme': 'Least square',
            'detection_method': 'swirling',
            'theoretical_model': 'Lamb-Oseen',
            'detection_threshold': 20,
            'correlation_threshold': 0.5,
            'boxsize': 6,
            'first': 5,
            'last': 6,
            'step': 1,
            'rmax': 0,
            'plot_method': 'fit',
            'xy_location': [0, 0],
            'num_cores': 1,
        }

        self.vfield = None
        self.vortices = None

        self.figures_list = []
        self.canvas_list = []

        self.current_vortex = 1
        self.total_vortices = None

        # Panneau principal
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)

        # Section des contrôles à droite (1/4 de la largeur)
        controls_panel = wx.Panel(self.main_tab)

        controls_sizer = wx.FlexGridSizer(rows=0, cols=2, vgap=5, hgap=10)
        controls_sizer.AddGrowableCol(1)  # Rendre la deuxième colonne extensible

        empty_label = wx.StaticText(controls_panel, label="")

        input_label = wx.StaticText(controls_panel, label="Input File:", name="input_file_label")
        self.input_text = wx.TextCtrl(controls_panel, value=self.analysis_parameters['input_filename'])
        self.file_status_icon = wx.StaticBitmap(controls_panel, bitmap=wx.BitmapBundle.FromBitmap(wx.NullBitmap))
        self.input_text.Bind(wx.EVT_TEXT, self.update_file_status)

        output_label = wx.StaticText(controls_panel, label="Output Directory:")
        self.output_text = wx.TextCtrl(controls_panel, value=self.analysis_parameters['output_directory'])

        meanfile_label = wx.StaticText(controls_panel, label="Mean file:", name="meanfile_label")
        self.meanfile_text = wx.TextCtrl(controls_panel, value=self.analysis_parameters['meanfile'])
        self.file_status_icon = wx.StaticBitmap(controls_panel, bitmap=wx.BitmapBundle.FromBitmap(wx.NullBitmap))
        self.meanfile_text.Bind(wx.EVT_TEXT, self.update_file_status)

        scheme_label = wx.StaticText(controls_panel, label="Scheme:")
        self.scheme_choice = wx.Choice(controls_panel, choices=["Second-order", "Least-square", "Fourth-order"])
        match self.analysis_parameters['scheme']:
            case 'Second-order':
                self.scheme_choice.SetSelection(0)
            case 'Fourth-order':
                self.scheme_choice.SetSelection(2)
            case _:
                self.scheme_choice.SetSelection(1)

        detection_label = wx.StaticText(controls_panel, label="Detection Method:")
        self.detection_choice = wx.Choice(controls_panel, choices=["Q", "delta", "swirling"])
        match self.analysis_parameters['detection_method']:
            case 'Q':
                self.detection_choice.SetSelection(0)
            case 'delta':
                self.detection_choice.SetSelection(1)
            case _:
                self.detection_choice.SetSelection(2)

        theoretical_model_label = wx.StaticText(controls_panel, label="Theoretical model:")
        self.theoretical_model_choice = wx.Choice(controls_panel, choices=["Rankine", "Lamb-Oseen", "Batchelor"])
        match self.analysis_parameters['theoretical_model']:
            case 'Rankine':
                self.theoretical_model_choice.SetSelection(0)
            case 'Batchelor':
                self.theoretical_model_choice.SetSelection(2)
            case _:
                self.theoretical_model_choice.SetSelection(1)

        file_type_label = wx.StaticText(controls_panel, label="File type:")


        self.file_type_choice = wx.Choice(
            controls_panel, choices=["PIV Tecplot (.dat)", "PIV NetCDF (.nc)", "DNS (.nc)", "OpenFOAM (.raw)", "HDF5 (.h5)"]
        )
        match self.analysis_parameters['file_type']:
            case 'PIV Tecplot (.dat)':
                self.file_type_choice.SetSelection(0)
            case 'PIV NetCDF (.nc)':
                self.file_type_choice.SetSelection(1)
            case 'DNS (.nc)':
                self.file_type_choice.SetSelection(2)
            case 'OpenFOAM (.raw)':
                self.file_type_choice.SetSelection(3)
            case 'HDF5 (.h5)':
                self.file_type_choice.SetSelection(4)

        max_radius_label = wx.StaticText(controls_panel, label="Max. radius:")
        self.max_radius_text = wx.TextCtrl(controls_panel, value=str(self.analysis_parameters['rmax']))

        threshold_label = wx.StaticText(controls_panel, label="Detection Threshold:")
        self.threshold_text = wx.TextCtrl(controls_panel, value=str(self.analysis_parameters['detection_threshold']))

        correlation_threshold_label = wx.StaticText(controls_panel, label="Correlation Threshold:")
        self.correlation_threshold_text = wx.TextCtrl(
            controls_panel, value=str(self.analysis_parameters['correlation_threshold'])
        )

        boxsize_label = wx.StaticText(controls_panel, label="Box Size:")
        self.boxsize_text = wx.TextCtrl(controls_panel, value=str(self.analysis_parameters['boxsize']))

        plot_method_label = wx.StaticText(controls_panel, label="Plot Method:")
        self.plot_method_choice = wx.Choice(controls_panel, choices=["fit", "detect", "fields", "xy location"])
        match self.analysis_parameters['plot_method']:
            case 'fit':
                self.plot_method_choice.SetSelection(0)
            case 'method':
                self.plot_method_choice.SetSelection(1)
            case 'fields':
                self.plot_method_choice.SetSelection(2)
            case 'xy location':
                self.plot_method_choice.SetSelection(3)

        xy_location_label = wx.StaticText(controls_panel, label="XY coordinates:")
        self.x_location = wx.TextCtrl(controls_panel, value=str(0))
        self.y_location = wx.TextCtrl(controls_panel, value=str(0))
        self.xy_location = [int(self.x_location.GetValue()), int(self.y_location.GetValue())]

        cores_label = wx.StaticText(controls_panel, label="Number of Cores:")
        max_cores = os.cpu_count()  # Get the maximum number of CPU cores
        self.cores_input_spin = wx.SpinCtrl(controls_panel, value='1', min=1, max=max_cores)

        self.start_button = wx.Button(controls_panel, label=" Start Analysis ")
        self.start_button.Bind(wx.EVT_BUTTON, self.on_start)

        # Ajout des éléments de contrôle dans la colonne à droite

        controls_sizer.Add(input_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.input_text, 0, wx.ALL | wx.EXPAND, 5)

        # Toggle Button for Time Series
        self.time_series_toggle = wx.ToggleButton(controls_panel, label="Toggle Time series")
        self.time_series_toggle.SetValue(False)
        self.time_series_toggle.Bind(wx.EVT_TOGGLEBUTTON, self.on_toggle_time_series)

        # Input fields for first, last, and step
        first_label = wx.StaticText(controls_panel, label="First:")
        self.first_text = wx.TextCtrl(controls_panel, value=str(self.analysis_parameters['first']))

        last_label = wx.StaticText(controls_panel, label="Last:")
        self.last_text = wx.TextCtrl(controls_panel, value=str(self.analysis_parameters['last']))

        step_label = wx.StaticText(controls_panel, label="Step:")
        self.step_text = wx.TextCtrl(controls_panel, value=str(self.analysis_parameters['step']))

        self.first_text.Enable(False)
        self.last_text.Enable(False)
        self.step_text.Enable(False)

        controls_sizer.Add(file_type_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.file_type_choice, 0, wx.ALL, 5)

        controls_sizer.Add(self.time_series_toggle, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(empty_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(first_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.first_text, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(last_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.last_text, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(step_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.step_text, 0, wx.ALL | wx.EXPAND, 5)

        controls_sizer.Add(output_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.output_text, 0, wx.ALL | wx.EXPAND, 5)
        
        controls_sizer.Add(meanfile_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.meanfile_text, 0, wx.ALL | wx.EXPAND, 5)
        
        controls_sizer.Add(scheme_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.scheme_choice, 0, wx.ALL, 5)
        
        controls_sizer.Add(detection_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.detection_choice, 0, 0, wx.ALL, 5)
        
        controls_sizer.Add(theoretical_model_label, 0, 0, wx.ALL, 5)
        controls_sizer.Add(self.theoretical_model_choice, 0, 0, wx.ALL, 5)
        
        controls_sizer.Add(threshold_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.threshold_text, 0, wx.ALL, 5)
        
        controls_sizer.Add(correlation_threshold_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.correlation_threshold_text, 0, wx.ALL, 5)
        
        controls_sizer.Add(max_radius_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.max_radius_text, 0, wx.ALL, 5)
        
        controls_sizer.Add(boxsize_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.boxsize_text, 0, wx.ALL, 5)
        
        controls_sizer.Add(plot_method_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_sizer.Add(self.plot_method_choice, 0, wx.ALL, 5)

        controls_xy_sizer = wx.FlexGridSizer(rows=0, cols=3, vgap=5, hgap=10)
        controls_sizer.Add(xy_location_label, 0, wx.ALL | wx.EXPAND, 5)
        controls_xy_sizer.Add(self.x_location, 0, wx.ALL, 5)
        controls_xy_sizer.Add(self.y_location, 0, wx.ALL, 5)
        controls_sizer.Add(controls_xy_sizer, 0, wx.ALL, 5)

        #controls_cores_sizer = wx.FlexGridSizer(rows=0, cols=2, vgap=5, hgap=10)
        controls_sizer.Add(cores_label, 0, wx.ALL, 5)
        controls_sizer.Add(self.cores_input_spin, 0, wx.ALL, 5)
        # controls_cores_sizer.Add(wx.StaticText(controls_panel, label="AA"), 0, wx.EXPAND)
        #controls_sizer.Add(controls_cores_sizer, 0, wx.ALL, 5)

        self.read_button = wx.Button(controls_panel, label="Read analysis parameters")
        self.save_button = wx.Button(controls_panel, label="Save analysis parameters")

        # Associer les boutons aux méthodes de lecture et sauvegarde
        self.read_button.Bind(wx.EVT_BUTTON, self.on_read_parameters)
        self.save_button.Bind(wx.EVT_BUTTON, self.on_save_parameters)

        controls_sizer.Add(self.read_button, 0, wx.ALL, 5)
        controls_sizer.Add(self.save_button, 0, wx.ALL, 5)


        controls_sizer.Add(wx.StaticText(controls_panel, label=" "), 0, wx.ALL, 10)
        
        controls_sizer.Add(self.start_button, 0, wx.ALL, 10)

        controls_panel.SetSizer(controls_sizer)

        figures_panel = wx.Panel(self.main_tab)
        figures_sizer = wx.GridSizer(rows=1, cols=2, gap=wx.Size(10, 10))

        self.flip_button_list = []
        self.update_button_list = []
        self.flip_states = []  # Pour stocker l'état du bouton flip pour chaque figure
        self.fig_bounds = []

        for i in range(2):  # Boucle pour créer les deux figures
            # Création d'une figure et de son canvas associé
            figure = plt.Figure()
            canvas = FigureCanvas(figures_panel, -1, figure)
            _ = figure.add_subplot(111)

            # Ajout de la figure et des contrôles dans le GridSizer
            figure_sizer = wx.BoxSizer(wx.VERTICAL)
            figure_sizer.Add(canvas, 1, wx.EXPAND)

            # Zones de texte pour xmin, xmax, ymin, ymax
            size_text = 50
            # Créez les champs de texte pour les axes
            xmin_text = wx.TextCtrl(figures_panel, value="-10", size=wx.Size(size_text, -1))
            xmax_text = wx.TextCtrl(figures_panel, value="10", size=wx.Size(size_text, -1))
            ymin_text = wx.TextCtrl(figures_panel, value="-10", size=wx.Size(size_text, -1))
            ymax_text = wx.TextCtrl(figures_panel, value="10", size=wx.Size(size_text, -1))

            # Stockez les champs pour un accès dans les autres méthodes
            self.xmin = xmin_text
            self.xmax = xmax_text
            self.ymin = ymin_text
            self.ymax = ymax_text

            self.fig_bounds.append([xmin_text, xmax_text, ymin_text, ymax_text])
            # Refresh button
            refresh_button = wx.Button(figures_panel, label="Refresh")
            refresh_button.Bind(wx.EVT_BUTTON, partial(self.on_update_axes, i))
            # Flip axis button
            flip_button = wx.ToggleButton(figures_panel, label="Flip Axis")
            self.flip_button_list.append(flip_button)
            self.flip_states.append(False)  # initialize flip value at false
            flip_button.Bind(wx.EVT_TOGGLEBUTTON, partial(self.on_flip_axis, i))

            # Layout pour les boutons sous chaque figure
            axis_controls = wx.BoxSizer(wx.HORIZONTAL)
            # Ajouter les champs de texte et le bouton flip axis dans le sizer horizontal
            size_panel = 5
            axis_controls.Add(wx.StaticText(figures_panel, label="xmin:"), 0, wx.ALL, size_panel)
            axis_controls.Add(xmin_text, 0, wx.ALL, size_panel)
            axis_controls.Add(wx.StaticText(figures_panel, label="xmax:"), 0, wx.ALL, size_panel)
            axis_controls.Add(xmax_text, 0, wx.ALL, size_panel)
            axis_controls.Add(wx.StaticText(figures_panel, label="ymin:"), 0, wx.ALL, size_panel)
            axis_controls.Add(ymin_text, 0, wx.ALL, size_panel)
            axis_controls.Add(wx.StaticText(figures_panel, label="ymax:"), 0, wx.ALL, size_panel)
            axis_controls.Add(ymax_text, 0, wx.ALL, size_panel)
            axis_controls.Add(refresh_button, 0, wx.ALL, size_panel)
            axis_controls.Add(flip_button, 0, wx.ALL, size_panel)

            # Save figure button
            save_button = wx.Button(figures_panel, label="Save Figure")
            save_button.Bind(wx.EVT_BUTTON, partial(self.on_save_figure, i))

            # Ajout des contrôles d'axes et du bouton Save Figure dans le sizer vertical
            figure_controls_sizer = wx.BoxSizer(wx.VERTICAL)
            figure_controls_sizer.Add(axis_controls, 0, wx.ALIGN_CENTER)
            display_vortex_sizer = wx.BoxSizer(wx.HORIZONTAL)
            display_vortex_sizer.Add(save_button, 0, wx.ALIGN_CENTER | wx.TOP, 5)

            if i == 1:
                # navigation button for right figure
                previous_button = wx.Button(figures_panel, label="←")
                next_button = wx.Button(figures_panel, label="→")
                self.current_vortex_text = wx.TextCtrl(
                    figures_panel, value="1", style=wx.TE_CENTER | wx.TE_PROCESS_ENTER
                )
                self.total_vortex_text = wx.StaticText(figures_panel, label="/ 1")

                # sizer for navigation between vortex
                display_vortex_sizer_spacing = 3
                display_vortex_sizer.Add(previous_button, 0, wx.ALL, display_vortex_sizer_spacing)
                display_vortex_sizer.Add(self.current_vortex_text, 0, wx.ALL, display_vortex_sizer_spacing)
                display_vortex_sizer.Add(self.total_vortex_text, 0, wx.ALL, display_vortex_sizer_spacing)
                display_vortex_sizer.Add(next_button, 0, wx.ALL, display_vortex_sizer_spacing)

                previous_button.Bind(wx.EVT_BUTTON, self.on_previous_vortex)
                next_button.Bind(wx.EVT_BUTTON, self.on_next_vortex)
                self.current_vortex_text.Bind(wx.EVT_TEXT_ENTER, self.on_vortex_text_enter)

            figure_controls_sizer.Add(display_vortex_sizer, 0, wx.ALIGN_CENTER | wx.TOP, 5)

            # Ajout de la figure et des contrôles dans le sizer de la figure
            figure_sizer.Add(figure_controls_sizer, 0, wx.EXPAND)

            figures_sizer.Add(figure_sizer, 1, wx.EXPAND)

            # store figures and canvas in lists
            self.figures_list.append(figure)
            self.canvas_list.append(canvas)

        figures_panel.SetSizer(figures_sizer)

        # Ajout des deux panels (contrôles et figures) au layout principal
        main_sizer.Add(figures_panel, 1, wx.EXPAND, 10)
        main_sizer.Add(controls_panel, 1, wx.EXPAND, 10)

        self.main_tab.SetSizer(main_sizer)

        self.SetTitle('Vortex Fitting - GUI')
        self.Centre()

        self.update_file_status()

    def on_toggle_time_series(self, _event):
        """
        Enable / disable the time series buttons

        :param _event: event associated with the button pressed.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        is_enabled = self.time_series_toggle.GetValue()
        # Activer ou désactiver les champs en fonction de l'état du toggle button
        self.first_text.Enable(is_enabled)
        self.last_text.Enable(is_enabled)
        self.step_text.Enable(is_enabled)
        self.update_file_status()

    def on_start(self, _event):
        """
        Get all the parameters and run the analysis

        :param _event: event associated with the button pressed.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        self.analysis_parameters['input_filename'] = self.input_text.GetValue()
        self.analysis_parameters['output_directory'] = self.output_text.GetValue()
        self.analysis_parameters['file_type'] = self.file_type_choice.GetStringSelection()
        self.analysis_parameters['meanfile'] = self.meanfile_text.GetValue()

        self.analysis_parameters['scheme'] = self.scheme_choice.GetStringSelection()
        self.analysis_parameters['detection_method'] = self.detection_choice.GetStringSelection()
        self.analysis_parameters['theoretical_model'] = self.theoretical_model_choice.GetStringSelection().lower()
        self.analysis_parameters['detection_threshold'] = float(self.threshold_text.GetValue())
        self.analysis_parameters['boxsize'] = int(self.boxsize_text.GetValue())
        self.analysis_parameters['rmax'] = float(self.max_radius_text.GetValue())

        self.analysis_parameters['first'] = int(self.first_text.GetValue())
        self.analysis_parameters['last'] = int(self.last_text.GetValue())
        self.analysis_parameters['step'] = int(self.step_text.GetValue())

        self.analysis_parameters['xy_location'] = [int(self.x_location.GetValue()), int(self.y_location.GetValue())]
        self.analysis_parameters['num_cores'] = self.cores_input_spin.GetValue()

        # flip_axis = self.flip_toggle_button.GetValue()  # Récupère l'état True/False du ToggleButton
        self.analysis_parameters['plot_method'] = self.plot_method_choice.GetStringSelection()

        print(f"Selected scheme: {self.scheme_choice.GetString(self.scheme_choice.GetSelection())}")
        print(f"Theoretical model: {self.analysis_parameters['theoretical_model']}")

        self.run_analysis(self.analysis_parameters)

    def on_update_axes(self, index, _event):
        """
        Method called when the "Refresh" button is pressed, to update the figure axes

        :param index: The index of the corresponding figure and canvas.
        :param _event: event associated with the button pressed.
        :type index: int
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        try:
            # Récupérez les valeurs des champs de texte
            xmin = float(self.fig_bounds[index][0].GetValue())
            ymin = float(self.fig_bounds[index][2].GetValue())
            xmax = float(self.fig_bounds[index][1].GetValue())
            ymax = float(self.fig_bounds[index][3].GetValue())

            # Mettre à jour les limites des axes
            self.figures_list[index].gca().set_xlim([xmin, xmax])
            self.figures_list[index].gca().set_ylim([ymin, ymax])

            # Redessiner la figure
            self.canvas_list[index].draw()
            self.refresh_canvas(self.canvas_list[index])
        except ValueError:
            # Ignorer les erreurs si les valeurs ne sont pas valides (par exemple, lors de la saisie)
            pass

    def on_flip_axis(self, index, _event):
        """
        Method called when the Flip Axis button is toggled for a specific figure.

        :param index: The index of the corresponding figure and canvas.
        :param _event: event associated with the button pressed.
        :type index: int
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        self.flip_states[index] = self.flip_button_list[index].GetValue()
        if index == 0:
            gui_plot_accepted(
                self.figures_list[index].gca(),
                self.canvas_list[index],
                self.vfield,
                self.vortices,
                self.vfield.detection_field,
                self.flip_states[index],
            )
        elif index == 1:
            gui_plot_vortex(
                self.figures_list[index].gca(),
                self.canvas_list[index],
                self.vfield,
                self.vortices[self.current_vortex - 1],
                self.analysis_parameters['theoretical_model'],
                self.flip_states[index],
            )

    def on_save_figure(self, index, _event):
        """
        Method called to save the current figure as a PNG file.

        :param index: The index of the figure to save.
        :param _event: event associated with the save button.
        :type index: int
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        # Open a file dialog to specify the save location
        with wx.FileDialog(
            self, "Save PNG file", wildcard="PNG files (*.png)|*.png", style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
        ) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  # User canceled the dialog

            # Save the figure to the specified path
            path = fileDialog.GetPath()
            self.figures_list[index].savefig(path, format="png")

    def on_previous_vortex(self, _event):
        """
        Method called when the "Previous" button is pressed to navigate to the previous vortex.

        :param _event: event associated with the button press.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        if self.current_vortex > 1:
            self.current_vortex -= 1
            self.current_vortex_text.SetValue(str(self.current_vortex))
            gui_plot_vortex(
                self.figures_list[1].gca(),
                self.canvas_list[1],
                self.vfield,
                self.vortices[self.current_vortex - 1],
                self.analysis_parameters['theoretical_model'],
                self.flip_states[1],
            )
            self.refresh_canvas(self.canvas_list[1])

    def on_next_vortex(self, _event):
        """
        Method called when the "Next" button is pressed to navigate to the next vortex.

        :param _event: event associated with the button press.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        if self.current_vortex < self.total_vortices:
            self.current_vortex += 1
            self.current_vortex_text.SetValue(str(self.current_vortex))
            gui_plot_vortex(
                self.figures_list[1].gca(),
                self.canvas_list[1],
                self.vfield,
                self.vortices[self.current_vortex - 1],
                self.analysis_parameters['theoretical_model'],
                self.flip_states[1],
            )
            self.refresh_canvas(self.canvas_list[1])

    def refresh_canvas(self, canvas):
        """
        Refresh the given canvas and its parent window.

        :param canvas: the matplotlib canvas to refresh.
        :type canvas: canvas

        :returns: empty
        :rtype: empty
        """
        canvas.draw()
        canvas.GetParent().Refresh()
        canvas.GetParent().Update()

    def on_vortex_text_enter(self, _event):
        """
        Method called when a new vortex number is entered in the text field.
        Validates the input and updates the vortex display accordingly.

        :param _event: event associated with pressing Enter in the text field.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        try:
            # Get the new vortex number from the text field
            new_vortex_number = int(self.current_vortex_text.GetValue())

            # Validate that the vortex number is within the allowed range
            if 1 <= new_vortex_number <= self.total_vortices:
                self.current_vortex = new_vortex_number
                gui_plot_vortex(
                    self.figures_list[1].gca(),
                    self.canvas_list[1],
                    self.vfield,
                    self.vortices[self.current_vortex - 1],
                    self.analysis_parameters['theoretical_model'],
                    self.flip_states[1],
                )
            else:
                wx.MessageBox(
                    f"Please set a number between 1 and {self.total_vortices}.", "Error", wx.OK | wx.ICON_ERROR
                )
        except ValueError:
            wx.MessageBox("Please set a correct number.", "Error", wx.OK | wx.ICON_ERROR)

    def update_file_status(self, event=None):
        """
        Check if the file exists and update the icon and input file label accordingly.

        :param event: event tu the update of the file status
        :type event: event

        :returns: empty
        :rtype: empty
        """

        file_path_std_file = self.input_text.GetValue()

        std_file_exists = os.path.isfile(file_path_std_file)
        time_file_exists = False
        if self.time_series_toggle.GetValue():
            time_file_exists = os.path.isfile(file_path_std_file.format(int(self.first_text.GetValue())))
        ok_std = std_file_exists or time_file_exists


        input_label = self.FindWindowByName("input_file_label")
        if input_label:
            input_label.SetLabel(f"Input file: {'✅' if ok_std else '❌'}")

        file_path_mean_file = self.meanfile_text.GetValue()
        mean_file_exists = os.path.isfile(file_path_mean_file)
        meanfile_label = self.FindWindowByName("meanfile_label")

        if meanfile_label:
            meanfile_label.SetLabel(f"Mean file: {'✅' if (mean_file_exists | (file_path_mean_file == '/')) else '❌'}")
        # Disable Start Analysis button if file not found
        self.start_button.Enable(ok_std)

        if event:
            event.Skip()

    def on_read_parameters(self, _event):
        """
        Method to read a configuration file and update the analysis parameters.

        :param _event: event associated with the file reading button.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        with wx.FileDialog(
            self,
            "Open configuration file",
            wildcard="Config files (*.cfg)|*.cfg",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST,
        ) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  # User canceled the dialog

            # Load the configuration file
            path = fileDialog.GetPath()
            config = configparser.ConfigParser()
            config.read(path)

            # Update the analysis parameters
            if 'Analysis' in config:
                self.analysis_parameters['scheme'] = config['Analysis'].get(
                    'scheme', self.analysis_parameters['scheme']
                )
                self.analysis_parameters['detection_method'] = config['Analysis'].get(
                    'detection_method', self.analysis_parameters['detection_method']
                )
                self.analysis_parameters['theoretical_model'] = config['Analysis'].get(
                    'theoretical_model', self.analysis_parameters['theoretical_model']
                )
                self.analysis_parameters['detection_threshold'] = float(
                    config['Analysis'].get('detection_threshold', self.analysis_parameters['detection_threshold'])
                )
                self.analysis_parameters['correlation_threshold'] = float(
                    config['Analysis'].get('correlation_threshold', self.analysis_parameters['correlation_threshold'])
                )
                self.analysis_parameters['boxsize'] = int(
                    config['Analysis'].get('boxsize', self.analysis_parameters['boxsize'])
                )
                self.analysis_parameters['first'] = int(
                    config['Analysis'].get('first', self.analysis_parameters['first'])
                )
                self.analysis_parameters['last'] = int(config['Analysis'].get('last', self.analysis_parameters['last']))
                self.analysis_parameters['step'] = int(config['Analysis'].get('step', self.analysis_parameters['step']))
                self.analysis_parameters['rmax'] = float(
                    config['Analysis'].get('rmax', self.analysis_parameters['rmax'])
                )
                self.analysis_parameters['input_filename'] = config['Analysis'].get(
                    'input_filename', self.analysis_parameters['input_filename']
                )
                self.analysis_parameters['output_directory'] = config['Analysis'].get(
                    'output_directory', self.analysis_parameters['output_directory']
                )
                self.analysis_parameters['meanfile'] = config['Analysis'].get(
                    'meanfile', self.analysis_parameters['meanfile']
                )
                self.analysis_parameters['plot_method'] = config['Analysis'].get(
                    'plot_method', self.analysis_parameters['plot_method']
                )

                # Ensure xy_location is correctly parsed
                xy_str = config['Analysis'].get('xy_location', "[0,0]")
                self.analysis_parameters['xy_location'] = [int(x) for x in xy_str.strip("[]").split(",")]
                self.analysis_parameters['num_cores'] = int(
                    config['Analysis'].get('num_cores', self.analysis_parameters['num_cores']))

            self.input_text.SetValue(self.analysis_parameters['input_filename'])
            self.output_text.SetValue(self.analysis_parameters['output_directory'])
            self.meanfile_text.SetValue(self.analysis_parameters['meanfile'])
            self.threshold_text.SetValue(str(self.analysis_parameters['detection_threshold']))
            self.correlation_threshold_text.SetValue(str(self.analysis_parameters['correlation_threshold']))
            self.boxsize_text.SetValue(str(self.analysis_parameters['boxsize']))
            self.first_text.SetValue(str(self.analysis_parameters['first']))
            self.last_text.SetValue(str(self.analysis_parameters['last']))
            self.step_text.SetValue(str(self.analysis_parameters['step']))
            self.max_radius_text.SetValue(str(self.analysis_parameters['rmax']))
            self.x_location.SetValue(str(self.analysis_parameters['xy_location'][0]))
            self.y_location.SetValue(str(self.analysis_parameters['xy_location'][1]))
            self.cores_input_spin.SetValue((str(self.analysis_parameters['num_cores'])))

            self.set_choice_selection(self.scheme_choice, self.analysis_parameters['scheme'])
            self.set_choice_selection(self.detection_choice, self.analysis_parameters['detection_method'])
            self.set_choice_selection(self.theoretical_model_choice, self.analysis_parameters['theoretical_model'])
            self.set_choice_selection(self.plot_method_choice, self.analysis_parameters['plot_method'])

            self.Update()
            self.Layout()
            self.Refresh()

            # print("Parameters loaded successfully:", self.analysis_parameters)

    def set_choice_selection(self, choice_widget, parameter_value):
        """
        Set choice selection based on a parameter value.

        :param choice_widget: the widget to update
        :type choice_widget: widget
        :param parameter_value: parameter from the Analysis structure
        :type parameter_value: string

        :returns: empty
        :rtype: empty
        """
        choices = [choice_widget.GetString(i) for i in range(choice_widget.GetCount())]

        for i, option in enumerate(choices):
            if option.lower() == str(parameter_value).lower():
                choice_widget.SetSelection(i)
                choice_widget.Update()
                choice_widget.Refresh()
                return
        # set the default values
        choice_widget.SetSelection(0)
        choice_widget.Update()
        choice_widget.Refresh()

    def on_save_parameters(self, _event):
        """
        Save the current analysis parameters to a configuration file (.cfg).

        :param _event: event associated with the save parameters button.
        :type _event: event

        :returns: empty
        :rtype: empty
        """
        # Update the analysis_parameters dictionary with the current UI values
        self.analysis_parameters['input_filename'] = self.input_text.GetValue()
        self.analysis_parameters['output_directory'] = self.output_text.GetValue()
        self.analysis_parameters['file_type'] = self.file_type_choice.GetStringSelection()
        self.analysis_parameters['meanfile'] = self.meanfile_text.GetValue()

        self.analysis_parameters['scheme'] = self.scheme_choice.GetStringSelection()
        self.analysis_parameters['detection_method'] = self.detection_choice.GetStringSelection()
        self.analysis_parameters['theoretical_model'] = self.theoretical_model_choice.GetStringSelection()
        self.analysis_parameters['detection_threshold'] = self.threshold_text.GetValue()
        self.analysis_parameters['correlation_threshold'] = self.correlation_threshold_text.GetValue()
        self.analysis_parameters['boxsize'] = self.boxsize_text.GetValue()
        self.analysis_parameters['first'] = self.first_text.GetValue()
        self.analysis_parameters['last'] = self.last_text.GetValue()
        self.analysis_parameters['step'] = self.step_text.GetValue()
        self.analysis_parameters['rmax'] = self.max_radius_text.GetValue()
        self.analysis_parameters['plot_method'] = self.plot_method_choice.GetStringSelection()
        self.analysis_parameters['xy_location'] = [int(self.x_location.GetValue()), int(self.y_location.GetValue())]
        self.analysis_parameters['num_cores'] = self.cores_input_spin.GetValue()

        with wx.FileDialog(
            self,
            "Save configuration file",
            wildcard="Config files (*.cfg)|*.cfg",
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        ) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  # L'utilisateur a annulé la boîte de dialogue

            # Récupérer le chemin sélectionné pour sauvegarder le fichier
            path = fileDialog.GetPath()

            # Créer l'objet ConfigParser et remplir la section Analysis
            config = configparser.ConfigParser()
            config['Analysis'] = {k: str(v) for k, v in self.analysis_parameters.items()}  # Ensure values are strings

            # Écriture des paramètres dans le fichier choisi
            with open(path, "w") as configfile:
                config.write(configfile)

            print(f"Parameters saved to {path}")

    def run_analysis(self, analysis_parameters):
        """
        Run the current analysis parameters.

        :param analysis_parameters: all the analysis parameters
        :type analysis_parameters: dict

        :returns: empty
        :rtype: empty
        """
        input_filename = analysis_parameters['input_filename']
        # output_directory = analysis_parameters['output_directory']
        file_type = analysis_parameters['file_type']
        meanfile = analysis_parameters['meanfile']
        scheme = analysis_parameters['scheme']
        detection_method = analysis_parameters['detection_method']
        theoretical_model = analysis_parameters['theoretical_model']
        detection_threshold = analysis_parameters['detection_threshold']
        correlation_threshold = analysis_parameters['correlation_threshold']
        boxsize = analysis_parameters['boxsize']
        first = analysis_parameters['first']
        last = analysis_parameters['last']
        step = analysis_parameters['step']
        rmax = analysis_parameters['rmax']
        plot_method = analysis_parameters['plot_method']
        xy_location = analysis_parameters['xy_location']
        num_cores = analysis_parameters['num_cores']

        for current_step in range(first, last + step, step):
            if self.time_series_toggle.GetValue():
                vfield = classes.VelocityField(input_filename, current_step, meanfile, file_type)
                output_time = current_step
            else:
                try:
                    vfield = classes.VelocityField(input_filename, 0, meanfile, file_type)

                except Exception as e:
                    wx.MessageBox(f"Error initializing VelocityField: {e} \n\n"
                                  f"Probably a wrong file type ?", "Error", wx.OK | wx.ICON_ERROR)
                    return None
                output_time = 0
            self.vfield = vfield

            # ---- DIFFERENCE APPROXIMATION ----#
            if scheme == "Fourth-order":
                vfield.derivative = schemes.fourth_order_diff(vfield)
            elif scheme == "Second-order":
                vfield.derivative = schemes.second_order_diff(vfield)
            elif scheme == "Least-square":
                vfield.derivative = schemes.least_square_diff(vfield)
            else:
                print('No scheme', scheme, 'found. Exiting!')
                sys.exit()
            # print(round(time.time() - lap,3), 'seconds')

            # ---- VORTICITY ----#

            vorticity = vfield.derivative['dvdx'] - vfield.derivative['dudy']

            detection_field = []
            if detection_method == 'Q':
                detection_field = detection.calc_q_criterion(vfield)
            elif detection_method == 'swirling':
                detection_field = detection.calc_swirling(vfield)
            elif detection_method == 'delta':
                detection_field = detection.calc_delta_criterion(vfield)
            vfield.detection_field = detection_field
            # print(round(time.time() - lap,3), 'seconds')

            if vfield.normalization_flag:
                print('Normalization for ', vfield.normalization_direction, ' direction')
                detection_field = fitting.normalize(detection_field, vfield.normalization_direction)  # normalization

            peaks = fitting.find_peaks(detection_field, detection_threshold, boxsize)
            vortices_counterclockwise, vortices_clockwise = fitting.direction_rotation(vorticity, peaks)

            vortices = list()

            if plot_method == 'xy location':
                x_location = int(xy_location[0])
                y_location = int(xy_location[1])
                detection_field = np.asarray(detection_field)
                detection_field_window = detection_field[
                    y_location - 10 : y_location + 10, x_location - 10 : x_location + 10
                ]
                x_index, y_index, u_data, v_data, _ = fitting.window(
                    vfield, x_location, y_location, 10, theoretical_model
                )
                fitting.plot_quiver(x_index, y_index, u_data, v_data, detection_field_window)

            if plot_method == 'fit':
                plt.close('all')
                vortices = fitting.get_vortices(
                    vfield, peaks, vorticity, rmax, correlation_threshold, theoretical_model, num_cores
                )
                print('---- Accepted vortices ----')
                print(len(vortices))
                self.total_vortices = len(vortices)
                if self.total_vortices < 1:
                    print("No vortex")
                    self.figures_list[1].gca().clear()
                    self.canvas_list[1].draw()
                else:
                    self.vortices = vortices
                    self.total_vortex_text.SetLabel('/ {:d}'.format(self.total_vortices))
                    args = argparse.Namespace(
                        detection_method=self.analysis_parameters['detection_method'],
                        scheme=self.analysis_parameters['scheme'],
                        box_size=self.analysis_parameters['boxsize'],
                        detection_threshold=self.analysis_parameters['detection_threshold'],
                        rmax=self.analysis_parameters['rmax'],
                        correlation_threshold=self.analysis_parameters['correlation_threshold'],
                        theoretical_model=self.analysis_parameters['theoretical_model'],
                        mean_filename=self.analysis_parameters['meanfile'],
                        file_type=self.analysis_parameters['file_type'],
                    )

                    output.create(self.analysis_parameters['output_directory'], args)
                    gui_plot_accepted(
                        self.figures_list[0].gca(),
                        self.canvas_list[0],
                        vfield,
                        vortices,
                        detection_field,
                        self.flip_states[0],
                    )
                    gui_plot_vortex(
                        self.figures_list[1].gca(),
                        self.canvas_list[1],
                        vfield,
                        vortices[self.current_vortex - 1],
                        theoretical_model,
                        self.flip_states[1],
                    )
                    print('Writing in: ', self.analysis_parameters['output_directory'])
                    output.write(vortices, self.analysis_parameters['output_directory'], output_time)

            if plot_method == 'detect':
                fitting.plot_detect(vortices_counterclockwise, vortices_clockwise, detection_field, False)
            if plot_method == 'fields':
                fitting.plot_fields(vfield, vorticity)
            vortex_info = "\n".join(
                [
                    f"{i+1}: r: {vortex[0]:.3f}, gamma: {vortex[1]:.2f}, xc: {vortex[2]:.2f}, yc: {vortex[3]:.2f}, correlation: {vortex[4]:.2f}, utheta: {vortex[8]:.2f}, uz0: {vortex[9]:.2f}"
                    for i, vortex in enumerate(vortices[:10])
                ]
            )
            if len(vortices) > 10:
                vortex_info += "\n... more ..."

            if not self.time_series_toggle.GetValue():
                wx.MessageBox(
                    f"Analysis completed. {self.total_vortices} vortices found.\n{vortex_info}",
                    "Success",
                    wx.OK | wx.ICON_INFORMATION,
                )

        return None


class UtilitiesTab(wx.Panel):
    """
    A tab in the GUI providing utilities for generating and converting NetCDF and ASCII files.

    Features:
    - Generate a NetCDF file with customizable parameters.
    - Convert a NetCDF file to an ASCII file.
    - Convert an ASCII file to a NetCDF file.
    """

    def __init__(self, parent):
        """
        Initializes the UtilitiesTab with GUI components for file generation and conversion.

        :param parent: parent window or panel
        :type parent: wx.Window

        :returns: empty
        :rtype: empty
        """
        super().__init__(parent)

        # Layout
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        # Title
        main_sizer.Add(wx.StaticText(self, label="Generate a netCDF example file"), 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Input fields
        self.inputs = {}
        fields = [
            ("core_radius", 5.0),
            ("gamma", 30.0),
            ("x_center", 64.0),
            ("y_center", 192.0),
            ("u_advection", 0.0),
            ("v_advection", 0.0),
            ("uz0", 0.5),
            ("ndim", 256),
            ("outfile", "test_netCDF.nc"),
            ("theoretical_model", "lamb-oseen"),
        ]
        dict = {
            'core_radius': "Core radius",
            'gamma': "Circulation Gamma",
            'x_center': "X center",
            'y_center': "Y center",
            'u_advection': "Advection velocity (u)",
            'v_advection': "Advection velocity (v)",
            'uz0': "Vortex center vertical velocity",
            'ndim': "Dimensions",
            'outfile': "Output file",
            'theoretical_model': "Theoretical model",
        }
        for field, default in fields:
            row_sizer = wx.BoxSizer(wx.HORIZONTAL)
            label = wx.StaticText(self, label=f"{dict[field]}: ")
            input_field = wx.TextCtrl(self, value=str(default))
            self.inputs[field] = input_field
            row_sizer.Add(label, 1, wx.EXPAND | wx.ALL, 5)
            row_sizer.Add(input_field, 2, wx.EXPAND | wx.ALL, 5)
            main_sizer.Add(row_sizer, 0, wx.EXPAND)

        # Generate button
        generate_button = wx.Button(self, label="Generate")
        generate_button.Bind(wx.EVT_BUTTON, self.on_generate)
        main_sizer.Add(generate_button, 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # Title for File Conversion
        main_sizer.Add(wx.StaticText(self, label="Convert files"), 0, wx.ALL | wx.ALIGN_CENTER, 10)

        # NetCDF to ASCII
        main_sizer.Add(wx.StaticText(self, label="NetCDF to ASCII"), 0, wx.ALL, 5)
        netcdf_to_ascii_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.netcdf_input = wx.TextCtrl(self, value="./test_netCDF.nc")
        self.ascii_output = wx.TextCtrl(self, value="./test_netCDF_cv.txt")
        browse_netcdf_btn = wx.Button(self, label="Browse")
        browse_netcdf_btn.Bind(wx.EVT_BUTTON, self.on_browse_netcdf)
        netcdf_to_ascii_sizer.Add(self.netcdf_input, 2, wx.EXPAND | wx.ALL, 5)
        netcdf_to_ascii_sizer.Add(self.ascii_output, 2, wx.EXPAND | wx.ALL, 5)
        netcdf_to_ascii_sizer.Add(browse_netcdf_btn, 1, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(netcdf_to_ascii_sizer, 0, wx.EXPAND)
        convert_netcdf_btn = wx.Button(self, label="Convert")
        convert_netcdf_btn.Bind(wx.EVT_BUTTON, self.convert_netcdf_to_ascii)
        main_sizer.Add(convert_netcdf_btn, 0, wx.ALL | wx.ALIGN_CENTER, 5)

        # ASCII to NetCDF
        main_sizer.Add(wx.StaticText(self, label="ASCII to NetCDF"), 0, wx.ALL, 5)
        ascii_to_netcdf_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ascii_input = wx.TextCtrl(self, value="./test_netCDF_cv.txt")
        self.netcdf_output = wx.TextCtrl(self, value="./test_netCDF_cv.nc")
        browse_ascii_btn = wx.Button(self, label="Browse")
        browse_ascii_btn.Bind(wx.EVT_BUTTON, self.on_browse_ascii)
        ascii_to_netcdf_sizer.Add(self.ascii_input, 2, wx.EXPAND | wx.ALL, 5)
        ascii_to_netcdf_sizer.Add(self.netcdf_output, 2, wx.EXPAND | wx.ALL, 5)
        ascii_to_netcdf_sizer.Add(browse_ascii_btn, 1, wx.EXPAND | wx.ALL, 5)
        main_sizer.Add(ascii_to_netcdf_sizer, 0, wx.EXPAND)
        convert_ascii_btn = wx.Button(self, label="Convert")
        convert_ascii_btn.Bind(wx.EVT_BUTTON, self.convert_ascii_to_netcdf)
        main_sizer.Add(convert_ascii_btn, 0, wx.ALL | wx.ALIGN_CENTER, 5)

        self.SetSizer(main_sizer)

    def on_generate(self, _event):
        """
        Generates a NetCDF file with parameters specified by the user in the input fields.

        The parameters include physical properties like core radius, circulation, and
        advection velocities, as well as output file name and grid dimensions.

        :param _event: event triggered when the "Generate" button is clicked.
        :type _event: wx._event

        :returns: empty
        :rtype: empty
        """
        params = {}
        for key, ctrl in self.inputs.items():
            value = ctrl.GetValue()
            match key:
                case "outfile":
                    params[key] = str(value)
                case "theoretical_model":
                    params[key] = str(value)
                case "ndim":
                    params[key] = int(value)
                case _:
                    params[key] = float(value)
        try:
            generateNetCDF.generateNetCDF(**params)
            wx.MessageBox("Fichier NetCDF généré avec succès.", "Info", wx.OK | wx.ICON_INFORMATION)
        except Exception as e:
            wx.MessageBox(f"Erreur : {e}", "Erreur", wx.OK | wx.ICON_ERROR)

    def on_browse_netcdf(self, _event):
        """
        Opens a file dialog to select a NetCDF file for conversion.

        The selected file path is displayed in the NetCDF input field.

        :param _event: event triggered when the "Browse" button for NetCDF is clicked.
        :type _event: wx._event
        """
        with wx.FileDialog(self, "Select NetCDF file", wildcard="NetCDF files (*.nc)|*.nc", style=wx.FD_OPEN) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                self.netcdf_input.SetValue(dlg.GetPath())

    def on_browse_ascii(self, _event):
        """
        Opens a file dialog to select an ASCII file for conversion.

        The selected file path is displayed in the ASCII input field.

        :param _event: event triggered when the "Browse" button for ASCII is clicked.
        :type _event: wx._event
        """
        with wx.FileDialog(self, "Select ASCII file", wildcard="Text files (*.txt)|*.txt", style=wx.FD_OPEN) as dlg:
            if dlg.ShowModal() == wx.ID_OK:
                self.ascii_input.SetValue(dlg.GetPath())

    def convert_netcdf_to_ascii(self, _event):
        """
        Converts a NetCDF file to an ASCII file.

        The NetCDF file path is taken from the NetCDF input field, and the output
        ASCII file path is specified in the ASCII output field.

        :param _event: event triggered when the "Convert" button for NetCDF to ASCII is clicked.
        :type _event: wx._event

        :returns: empty
        :rtype: empty
        """
        netcdf_path = self.netcdf_input.GetValue()
        ascii_path = self.ascii_output.GetValue()
        wx.MessageBox(f"Converting {netcdf_path} to {ascii_path}", "Info", wx.OK)
        convertToASCII.netcdf_to_ascii(netcdf_path, ascii_path)

    def convert_ascii_to_netcdf(self, _event):
        """
        Converts an ASCII file to a NetCDF file.

        The ASCII file path is taken from the ASCII input field, and the output
        NetCDF file path is specified in the NetCDF output field.

        :param _event: event triggered when the "Convert" button for ASCII to NetCDF is clicked.
        :type _event: wx._event
        """
        ascii_path = self.ascii_input.GetValue()
        netcdf_path = self.netcdf_output.GetValue()
        wx.MessageBox(f"Converting {ascii_path} to {netcdf_path}", "Info", wx.OK)
        convertToNC.ascii_to_netcdf(ascii_path, netcdf_path)


class AboutTab(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)

        main_sizer = wx.BoxSizer(wx.VERTICAL)

        about_button = wx.Button(self, label="About VortexFitting")
        about_button.Bind(wx.EVT_BUTTON, self.show_about)

        main_sizer.AddStretchSpacer()
        main_sizer.Add(about_button, 0, wx.ALIGN_CENTER | wx.ALL, 10)
        main_sizer.AddStretchSpacer()

        self.SetSizer(main_sizer)

    def show_about(self, _event):
        """
        Define an "About" window, with information about VortexFitting.

        :param _event: event triggered when the "About" button, on About panel, is clicked.
        :type _event: wx._event

        :returns: empty
        :rtype: empty
        """
        info = wx.adv.AboutDialogInfo()

        info.SetName("VortexFitting - GUI")
        info.SetVersion("2.0.0")

        info.SetDescription(
            "VortexFitting is a fluid mechanics post-processing tool "
            "developed in Python.\n\n"
            "It detects vortices in a flow and evaluates their properties.\n\n"
            "Available models:\n"
            "- Rankine\n"
            "- Lamb-Oseen\n"
            "- Batchelor\n\n"
            "Detection methods:\n"
            "- Swirling strength\n"
            "- Delta criterion\n"
            "- Q criterion\n\n"
            "Numerical schemes:\n"
            "- Second-order\n"
            "- Fourth-order\n"
            "- Least-square filter\n\n"
            "Compatible formats:\n"
            "- netCDF\n"
            "- HDF5\n"
            "- raw/dat\n\n"
            "Parallel processing available."
        )

        info.SetWebSite(
            "https://github.com/guilindner/VortexFitting",
            "VortexFitting on GitHub"
        )


        info.AddDeveloper("Guilherme Linder")
        info.AddDeveloper("Yann Devaux")
        info.AddDeveloper("Ilkay Solak")

        info.SetLicence(
                "Original publication on SoftwareX:\n"
                "https://doi.org/10.1016/j.softx.2020.100604\n\n"
                "Creative Commons license"
        )

        wx.adv.AboutBox(info)


class MyApp(wx.App):
    def OnInit(self):
        frame = MainFrame(None)
        frame.Show()
        return True


if __name__ == "__main__":
    main()
