import json
import numpy as np
import shapely as sp
import gdspy
import matplotlib.pyplot as plt
import yaml as yml
import sys
import os
import klayout.db as db

# Read the polygons saved to file, after the initial Julia script


# Opening the yml file, as to read the coordinates and the name of the GDS file, that we wish to clip.
def cut_gds(coordinates, file):
    with open("DiaQ.yml", "r") as stream:
        try:
            data = yml.safe_load(stream)
        except yml.YAMLError as exc:
            print(exc)

    x_1 = coordinates['x1']
    y_1 = coordinates['y1']
    x_2 = coordinates['x2']
    y_2 = coordinates['y2']
    name = file['name']

    # Creating a new layout class

    ly = db.Layout()

    # Reading the GDS file

    ly.read(name + ".gds")

    # Designating the layer that belongs to the clip box

    clip_layer = ly.layer(100, 0)

    # Note that we use "bbox" to get the integer-unit bounding box

    clip_box = ly.top_cell().bbox_per_layer(clip_layer)

    # Creating the clip-box, which is designated by the coordinates that we read from the yml file

    b2 = db.DBox(x_1, y_1, x_2, y_2)

    # The clip method only takes integer-unit clip boxes so far.
    # It takes and returns cell indexes.

    clip_cell_index = ly.clip(ly.top_cell().cell_index(), b2)

    # Creating a new cell, which is the clipped cell and naming it.

    clip_cell = ly.cell(clip_cell_index)
    clip_cell.name = "SLICED"

    # Reading the bounding box of the clipped cell, such that we can transform it. 
    # This is done, so that no matter where in the design the clip box is placed, the clipped cell will be placed around the origin.
    # When the database-unit is small, such as for nanometer designs, the int32 type, which is the default data-type in klayout,
    # is perhaps not capable of representing the large numbers of coordinates,
    # especially when scaled up.

    bbox1 = clip_cell.bbox()

    # Defining a transformation via the Trans class.

    trans = db.Trans(-bbox1.p1.x + 200, -bbox1.p1.y + 200)

    # Actually transforming the cell.

    clip_cell.transform(trans)

    # Creating a SaveLayoutOptions class, which is used to specify various things when saving,
    # such as the scale factor and the database-unit.

    sAve = db.SaveLayoutOptions()
    sAve.scale_factor = 1
    sAve.dbu = 0.001

    # Saving the clipped cell as a new GDS file.

    clip_cell.write(name + "_sliced.gds", sAve)

def read_yaml(file_name):
    
    """
    :file_name: The name of the file, which contains the gate names, with the .yml file extension. Do not input the file extension.

    Returns:
        _name_dict_: A dictionary of the gate names
    """
    
    yaml_file_name = file_name + ".yml"
    with open(yaml_file_name, "r") as stream:
        try:
            conf = yml.safe_load(stream)
        except yml.YAMLError as exc:
            print(exc)

    name_dict = conf["gate-names"]
    diagnostics_values = conf["diagnostics-values"]
    gds_coordinates = conf["gds-coordinates"]
    gds_file = conf["gds-file"]
    
    return name_dict, diagnostics_values, gds_coordinates, gds_file

def path_checking():
    
    """
    Checks if the Inter-language files folder exists. If it does not, it creates the folder.
    Reads the conf.yml file
    """

    path = os.path.isdir('Inter-language files')
        
    if path == False:
        cwd = os.getcwd()
        inter_file_path = os.path.join(cwd, 'Inter-language files')
        os.mkdir(inter_file_path)
        
    with open("conf.yml", "r") as stream:
        try:
            conf = yml.safe_load(stream)
        except yml.YAMLError as exc:
            print(exc)


def import_gds(file_name):
    
    """
    :file_name: The name of the file, which contains the design, with the .gds file extension
    Returns:
        :polygons_init: a list of the polygons read in the design file
        :layers: a list of the layers of the polygons 
    """

    read_file = gdspy.GdsLibrary(infile=file_name);

    main_cell = read_file.top_level()[0]

    full_dict = main_cell.get_polygons(by_spec = True)
    
    polygons_init = []

    layers = []

    for key in full_dict:
        for i in range(len(full_dict[key])):
            polygon_temp = sp.geometry.Polygon(full_dict[key][i])
            polygons_init.append(polygon_temp)
            layers.append(key[0])

    return polygons_init, layers

def gate_naming(polygons_init):
    
    """
    
    An interactive function, which allows the user to name the gates. The function will display the gates, and the user
    inputs the name of the gate in the terminal. The function will then save the gate names to a .json file, which can be
    read by the populate_gates function.
    
    :polygons_init: A list of the polygons, which are the gates
    :return: A dictionary of the gate names and the corresponding gate objects
    """

    gate_names = {}
    count = 0
    plt.ion() 
    fig, ax = plt.subplots()
    plt.show()
    for i in polygons_init:
        
        for j in polygons_init:    

            x, y = j.exterior.xy
            if i == j:
                ax.fill(x, y, color = 'red')
            else:
                ax.fill(x, y, color='#6699cc', alpha = 0.7)

        gate_name = input("Please input the name of the gate: \n")
        
        gate_names.update({f'{count}': gate_name})

        ax.clear()

        count += 1
        
    with open("Inter-language files/gate_names.json", "w") as fp:
        json.dump(gate_names, fp) 

def populate_gates(max_distance, naming = False):
    
    """
    :max_distance: The maximum distance between two gates, for them to be considered neighbours
    :return: A dictionary of the gate names and the corresponding gate objects
    """

    path_checking()

    polygons_init, layers = import_gds('AXL MkVI QT828_sliced.gds')

    if naming:

        print("Do you want to name the gates? (Y/n)")
        
        if input() == "Y":
            gate_naming(polygons_init)

    try:
        with open("Inter-language files/gate_names.json", "r") as read_file:
            gate_names = json.load(read_file)
    except:
        print("The file gate_names.json does not exist in the Inter-language files folder. Please run the gate_naming function.")
        sys.exit()

    gates = {v : Gate(v, layers[i], polygons_init[i]) for i, (k, v) in enumerate(gate_names.items())}

    for i_index in gate_names:
        i = gate_names[f'{i_index}']
        for j_index in gate_names:
            j = gate_names[f'{j_index}']
            inter = sp.intersection(gates[f'{i}'].polygon, gates[f'{j}'].polygon)
            dist = sp.distance(gates[f'{i}'].polygon, gates[f'{j}'].polygon)
            if inter.area > 0 and i != j:
                gates[f'{i}'].update_overlaps(gates[f'{j}'], inter.area)
            if dist < max_distance and i != j and gates[f'{i}'].layer == gates[f'{j}'].layer:
                gates[f'{i}'].update_neighbour(gates[f'{j}'], dist)
                
    return gates



class Leakage_Matrix:
    # Defining the inherent properties of the class. These will be saved internally in the class, and can only be
    # manipulated through the class functions. This makes it easier to make agent-based simulations.
    
    def __init__(self, leakage_matrix, threshold, name_dict):
        
        """
        :leakage_matrix: A matrix of the leakage values between the gates
        :threshold: The threshold for when the numerical value in the matrix is considered
        to be actually leaking
        :name_dict: A dictionary, which maps the gate names to the indices of the leakage matrix
        """
        
        self.leakage_matrix = leakage_matrix
        self.threshold = threshold
        self.name_dict = name_dict
        
    def to_binary(self):
        """
        :return: A binary matrix, where the elements are 1 if the leakage is above the threshold and 0 if it is below
        """
        
        return np.where(self.leakage_matrix > self.threshold, 1, 0)
    
    def failure_modes(self, gate_dict, show_plot):
        """
        :binary: A binary matrix, where the elements are 1 if the leakage is above the threshold and 0 if it is below
        :gate_dict: A dictionary, which maps the gate names to the gate objects
        :show_plot: A boolean, which determines if the plot of the failure modes should be shown
        :return: A histogram of the failure modes
        """
    
        ohmics = []
    
        #Remember to name the ohmics in gate_names.yml file with the key "ohmic_i", where i is the index of the ohmic
        # otherwise this algorithm will not work
    
        for key in self.name_dict:
            if 'ohmic' in key:
                ohmics.append(self.name_dict[key])
    
        binary = self.to_binary()
        
        index_array = np.matrix.nonzero(binary)

        failure_modes = []

        name_values = list(self.name_dict.values())
        
        for index in range(len(index_array[0])):
            x, y = index_array[0][index], index_array[1][index]
            gate_x = gate_dict[name_values[x]]
            gate_y = gate_dict[name_values[y]]
            x_neighbours = list(gate_x.neighbours.keys())
            y_neighbours = list(gate_y.neighbours.keys())
            if gate_x.name in y_neighbours and gate_y.name in x_neighbours:
                failure_modes.append("Lift-off")
                gate_x.update_leakage(gate_y, "Lift-off")
                gate_y.update_leakage(gate_x, "Lift-off")
            x_overlap = list(gate_x.overlap.keys())
            y_overlap = list(gate_y.overlap.keys())
            if gate_x.name in y_overlap and gate_y.name in x_overlap:
                if gate_x.name in ohmics or gate_y.name in ohmics:
                    failure_modes.append("ALD/Substrate")
                    gate_x.update_leakage(gate_y, "ALD/Substrate")
                    gate_y.update_leakage(gate_x, "ALD/Substrate")
                else:
                    failure_modes.append("ALD")
                    gate_x.update_leakage(gate_y, "ALD")
                    gate_y.update_leakage(gate_x, "ALD")
            if gate_x.name in ohmics or gate_y.name in ohmics:
                failure_modes.append("Substrate")
                gate_x.update_leakage(gate_y, "Substrate")
                gate_y.update_leakage(gate_x, "Substrate")
                
                
        if show_plot:      
        
            fig, axs = plt.subplots(1, 1, figsize=(10, 5), tight_layout=True)
            
            keys, counts = np.unique(failure_modes, return_counts=True)

            patches = axs.bar(keys, counts, facecolor = 'gray', edgecolor = 'red')
            
            
            
            bin_centers = [patches[i].get_x() + 0.5*patches[i].get_width() for i in range(0, len(patches))]
            counts = [patches[i].get_height() for i in range(0, len(patches))]
            
            height = 0.025*np.max(counts)
            
            for i in range(0, len(counts)):
                y_pos = counts[i] + height
                label = str(counts[i]) 
                axs.text(bin_centers[i] - 0.05, y_pos, label)
            
            axs.set_title("Failure modes")
            
            axs.set_ylabel("Frequency")
            
            plt.show()


class Gate:
    # Defining the inherent properties of the class. These will be saved internally in the class, and can only be
    # manipulated through the class functions. This makes it easier to make agent-based simulations.
    def __init__(self, name, layer, polygon):
        
        """
        :name: The name of the gate
        :layer: The layer of the gate
        :polygon: The polygon of the gate
        :overlap: A dictionary of the gates, which it overlaps with
        :neighbours: A dictionary of the gates, which are neighbours
        :leakage: A dictionary of the gates, which it is leaking to
        """
        
        self.name = name
        self.layer = layer
        self.polygon = polygon
        self.overlap = dict()
        self.neighbours = dict()
        self.leakage = dict()
        

    # Definings the various functions, which updates the internal properties of the class.


    def update_overlaps(self, other, area):
        """
        :other: the other gate it overlaps with
        :area: the area of the overlap
        :return: Does not return anything, but changes the dictionary of overlappings gates and adds the overlapping area
        """
        
        self.overlap.update({other.name: area})

    def update_neighbour(self, other, distance):
        """_summary_

        Args:
            other (Gate): The other gate, which is a neighbour
            distance (Length): The length between the two gates
        
        :return: Does not return anything, but changes the dictionary of neighbouring gates and adds the distance
        """
        
        self.neighbours.update({other.name: distance})
        
    def update_leakage(self, other, mode):
        """
        :other: the other gate, which it is leaking ti
        :mode: the mode of the leakage
        :return: Does not return anything, but changes the dictionary of leaking gates and adds the mode
        """
        
        self.leakage.update({other.name: mode})
        


# gates = populate_gates(max_distance = 100)

# leakage_matrix_init = np.random.random((29, 29))

# name_dict = {"0": "P6", "1": "P5", "2": "P4", "3": "SP2", "4": "SP1", "5": "P1", "6": "P2", "7": "P3", "8": "SP2PB", "9": "B7", "10": "B6", "11": "B5", "12": "B4", "13": "RTB", "14": "LTB", "15": "B3", "16": "RBB", "17": "SP1PB", "18": "B2", "19": "B1", "20": "LBB", "21": "S2", "22": "SS2", "23": "SS1", "24": "S1", "25": "RTO", "26": "LTO", "27": "RBO", "28": "LBO"}# "30": 5, "31": 6, "32": 1, "33": 2, "34": 3, "35": 4, "36": 5, "37": 6,
#              #"38": 1, "39": 2, "40": 3, "41": 4, "42": 5, "43": 6, "44": 1}

# leakage_matrix = Leakage_Matrix(leakage_matrix_init, 0.9, name_dict)

# leakage_matrix.failure_modes(gates)