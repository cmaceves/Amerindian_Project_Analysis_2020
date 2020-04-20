from sklearn.svm import LinearSVC 
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
import numpy as np
import pandas as pd 
import seaborn as sns
import math
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import matplotlib
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import argparse
from sklearn.decomposition import PCA
import seaborn as sns
import squarify
URBANIZATION_COLORS =  ["#0D0887","#9C179E","#ED7953","#F0F921", "#d2d2d2"]

###function takes in the annotation page and collapses based on the SMILES/INCHI/Name
def collapse_annotations(df):
    df.drop(["row retention time", "row m/z", "row ID"], axis=1, inplace=True)   
    df['scans'] = df["scans"].astype(str)
    df['scans'] = df['scans'].str.split(', ')

    #first we collapse the dataframe based on identical compound names
    pre_compound_df = df.drop(["SMILES","INCHI"], axis=1)
    columns = pre_compound_df.columns.tolist()
    sum_dict = {k:'sum' for k in columns}
    compound_collapsed = pre_compound_df.groupby('Compound_Name').agg(sum_dict)
    
    #collapse on smiles 
    full_compounds = df["Compound_Name"].tolist()
    full_smiles = df["SMILES"].tolist()
    
    new_row =[]

    for index, row in compound_collapsed.items():
        for compound_index, value in row.items():
            index_list = [i for i, x in enumerate(full_compounds) if x == compound_index]
            all_smiles_values = [full_smiles[i] for i in index_list if str(full_smiles[i]) != 'nan']
            if len(all_smiles_values) == 0:
                new_row.append('N/A')
            else:
                new_row.append(all_smiles_values[0])
        break
   

    compound_collapsed["SMILES"] = new_row
    compound_collapsed['SMILES'] = compound_collapsed['SMILES'].replace({'N/A': np.nan})
    compound_collapsed_df2 = compound_collapsed[compound_collapsed['SMILES'].notnull()]
    nan_compound_collapsed = compound_collapsed[compound_collapsed['SMILES'].isnull()]
    columns = compound_collapsed_df2.columns.tolist()
    sum_dict = {k:'sum' for k in columns}
    smiles_collapsed = compound_collapsed_df2.groupby('SMILES').agg(sum_dict)
    smiles_collapsed = pd.concat([smiles_collapsed,nan_compound_collapsed])
    smiles_collapsed.drop("SMILES", axis=1, inplace=True)

    #resetting the index to compound name
    smiles_collapsed_compounds = smiles_collapsed["Compound_Name"].tolist()
    scans = smiles_collapsed["scans"].tolist()
    new_index = [] 
    for scan,string in zip(scans, smiles_collapsed_compounds):
        divides = len(scan)
        firstpart = string[:len(string)//divides]
        new_index.append(firstpart)
    smiles_collapsed["Compound_Name"] = new_index
    smiles_collapsed.set_index("Compound_Name", inplace=True)
    smiles_collapsed["Compound_Name_2"] = new_index
    
    #collapse on inchi
    full_inchi = df["INCHI"].tolist()
    new_row =[]

    for index, row in smiles_collapsed.items():
        for compound_index, value in row.items():
            index_list = [i for i, x in enumerate(full_compounds) if x == compound_index]
            all_inchi_values = [full_inchi[i] for i in index_list if str(full_inchi[i]) != 'nan'] 
            if len(all_inchi_values) == 0:
                new_row.append('N/A')
            else:
                new_row.append(all_inchi_values[0])
        break
   
    smiles_collapsed["INCHI"] = new_row 
    smiles_collapsed['INCHI'] = smiles_collapsed['INCHI'].replace({'N/A': np.nan})
    smiles_collapsed_df2 = smiles_collapsed[smiles_collapsed['INCHI'].notnull()]
    nan_smiles_collapsed = smiles_collapsed[smiles_collapsed['INCHI'].isnull()]
    columns = smiles_collapsed_df2.columns.tolist()
    sum_dict = {k:'sum' for k in columns}
    inchi_collapsed = smiles_collapsed_df2.groupby('INCHI').agg(sum_dict)
    inchi_collapsed = pd.concat([inchi_collapsed, nan_smiles_collapsed]) 
    
    
    #resetting the index to compound name
    inchi_collapsed_compounds = inchi_collapsed["Compound_Name_2"].tolist()
    scans = inchi_collapsed["scans"].tolist()
    new_index = []
    #for string in inchi_collapsed_compounds:
    #    firstpart, secondpart = string[:len(string)//2], string[len(string)//2:]
    #    new_index.append(firstpart)
    #inchi_collapsed["Compound_Name_2"] = new_index
    inchi_collapsed.set_index("Compound_Name_2", inplace=True)
    print(inchi_collapsed)
    return(inchi_collapsed)


###function prints out color bar for color maps
def colorbar_colormap(num_items):
    viridis = cm.get_cmap('viridis', 256)
    
    viridisBig = cm.get_cmap('viridis', num_items)
    newcmp = ListedColormap(viridisBig(np.linspace(0.25, 0.75, 256)))

    colormaps = [viridis, newcmp]
    
    np.random.seed(19680801)
    data = np.random.randn(30, 30)
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            constrained_layout=True, squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=0, vmax=100)
        fig.colorbar(psm, ax=ax)
    plt.show()

###outputs boxplots to track univariate distributions
def pdf(matrix, labels, features, scan, compound_name=None): 
    cmap = matplotlib.cm.get_cmap(name="plasma")
    #cmap = ["red", "orange", "yellow", "lime", "green", "cyan", "blue", "purple", "lavender"]
    location = features.index(scan)
        
    for i,row in enumerate(matrix.T):
        if i == location :
            new_array = row
            break
    print(new_array)
    unique_labels = np.unique(labels).tolist()
    unique_counter = np.unique(labels, return_counts=True) 
    value_array = [0.0] * len(unique_counter[0])    

    count_occurance_1 = []
    count_occurance_2= []
    count_occurance_3 = []
    count_occurance_4 = []

    for item1, item2 in zip(new_array, labels):
        if item1 != 0:
            if item2 == 1:
                count_occurance_1.append(item1)
            if item2 == 2:
                count_occurance_2.append(item1)
            if item2 == 3:
                count_occurance_3.append(item1) 
            if item2 == 4:
                count_occurance_4.append(item1)

    #log scale the data
    inten_array = []
    for number in new_array:
        if number != 0:
            inten_array.append(math.log10(number) * -1)
        else:
            inten_array.append(0)

    print(len(labels))
    print(len(inten_array))

    #create a dataframe
    temp_df = pd.DataFrame({"labels":labels, "values":inten_array})

    #intialize figure and subplot
    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot()

    #create boxplot and swarmplot
    ax=sns.boxplot(y='values', x='labels', data=temp_df,hue=labels, palette= URBANIZATION_COLORS, dodge =False, fliersize=0)
    ax=sns.swarmplot(x=labels, y = inten_array, color ='black')
    ax.legend_.remove()
    #sns.boxplot(count_occurance_2, linewidth=1.5, label='Medium', color = URBANIZATION_COLORS[1])
    #sns.lineplot(count_occurance_3, linewidth=1.5, label='High', color = URBANIZATION_COLORS[2])
    #sns.lineplot(count_occurance_4, linewidth=1.5, label='Mestizo', color = URBANIZATION_COLORS[3])
    
    ax.set_ylabel('Intensity', fontsize=15)
    ax.set_xlabel('Group', fontsize=15)
    ax.set_xticklabels(['Low', "Medium", 'High', 'Mestizo'],fontsize=15)
    ax.tick_params(labelsize=15)
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='darkgray')
    ax.xaxis.grid(color='darkgray')    
    ax.figure.savefig("./boxplot_urine_%s.pdf" %compound_name)


#meant to calculate the total sum of squares
def calcaulte_SST(x_coordinates, y_coordinates):
    x_average = np.mean(x_coordinates)
    y_average = np.mean(y_coordinates)
    
    ss_total = 0
    
    for x,y in zip(x_coordinates, y_coordinates):
        ss_x = (x - x_average) ** 2
        ss_y = (y - y_average) ** 2

        ss_total += ss_x
        ss_total += ss_y
    
    return(ss_total, x_average, y_average)

#meant to calculate the sum of squares 
def calculate_SSW(cluster_x, cluster_y):
    x_average = np.mean(cluster_x)
    y_average = np.mean(cluster_y)
    
    ss_total = 0

    for x,y in zip(cluster_x, cluster_y):
        ss_x = (x - x_average) ** 2 
        ss_y = (y - y_average) ** 2
                    
        ss_total += ss_x
        ss_total += ss_y
 
    return(ss_total)

def calculate_RMSSTD(cluster_x, cluster_y):
    x_average = np.mean(cluster_x)
    y_average = np.mean(cluster_y)

    ss_x = 0
    ss_y = 0

    dof_x = len(cluster_x) - 1
    dof_y = len(cluster_y) - 1

    for x,y in zip(cluster_x, cluster_y):
        ss_x += (x - x_average) ** 2
        ss_y += (y - y_average) ** 2

    RMSSTD = (ss_x + ss_y) / (dof_x + dof_y)
    RMSSTD = math.sqrt (RMSSTD)

    return(RMSSTD)

###function row sum normalizes each samples
def normalize_quantitative(matrix):
    total = np.sum(matrix, axis = 0)
    matrix = np.divide(matrix, total)
    matrix[matrix<=0.00005]=0
    matrix[matrix>0.00005]=1
    return(matrix)

###defines the function where we split the data into a 70/30 scheme for training/prediction respectively
### shuffle data
def split_data(df, metadata, metadata_category, body_site = None):
    #filtering out things that don't have the metadata we want
    #etadata = metadata[(metadata[metadata_category] != "not specified")]     
    
    if body_site != None:
        metadata = metadata[(metadata["Body_site"] == body_site)]
        
    metadata.set_index("Name_MS", inplace=True) 
   
    droppy_cols = metadata.columns.tolist()
    droppy_cols = [item for item in droppy_cols if item.startswith("Unnamed")]
    metadata.drop(droppy_cols, axis=1, inplace=True)

    metadata = metadata.T

    droppy_cols = metadata.columns.tolist()
    droppy_cols = [item for item in droppy_cols if item.startswith("Unnamed")]
    metadata.drop(droppy_cols, axis=1, inplace=True)

    df, metadata = df.align(metadata, join = "inner", axis =1)
    metadata = metadata.T
    
    #mapping condition for a specific metadata category
    if metadata_category == "Global_individual_urbanization_groups":
        groups={'Low.ind': 1, 'Medium.ind': 2, 'High.ind': 3, 'High.ind_mes': 4, "not specified" : "nan", "not applicable" : "nan"}    
        metadata["Global_individual_urbanization_groups"]=[groups[i] for i in metadata["Global_individual_urbanization_groups"].tolist()]
    
    metaSamples=metadata.index.tolist()
    globalUrb = metadata["%s" %metadata_category].tolist()
    samples = df.columns.tolist()

    zippy = [(item1, item2) for item1,item2 in zip(metaSamples,globalUrb) if item1 in samples]
    masterList = [item[0] for item in zippy]
    globalUrb = [item[1] for item in zippy]

    print("samples", len(samples))
    print("metaSamples", len(metaSamples))
    print("master", len(masterList))

    matrix = df.to_numpy() 
    
    indices_to_drop = [count for count,item in enumerate(samples) if item not in masterList]
    matrix=np.delete(matrix,indices_to_drop, 1)
    labels=np.array([globalUrb])
       
    return(np.transpose(matrix), labels)

###function generates a bar plot    
def bar_chart_general(matrix, labels, features, scan, compound_name = None):
    rows_to_collapse = []
    new_array = [0.0]*matrix.shape[0]
    
    cmap = matplotlib.cm.get_cmap(name="plasma")
    #cmap = ["red", "orange", "yellow", "lime", "green", "cyan", "blue", "purple", "lavender"]
    location = features.index(scan)
        
    for i,row in enumerate(matrix.T):
        if i == location :
            new_array = row
            break
    
    unique_labels = np.unique(labels).tolist()
    unique_counter = np.unique(labels, return_counts=True) 
    value_array = [0.0] * len(unique_counter[0])

    count_occurance_1 = 0
    count_occurance_2= 0
    count_occurance_3 = 0
    count_occurance_4 = 0

    for item1, item2 in zip(new_array, labels):
        if item1 != 0:
            locate_position = unique_labels.index(item2)
            value_array[locate_position] += 1 
    
    list_percent = [x/y for x,y in zip(value_array,unique_counter[1])]
    
    fig = plt.figure()
    ax=fig.add_subplot()
    
    ax.bar(unique_labels, list_percent, color=URBANIZATION_COLORS) 
    #plt.xticks(rotation=45)
    #ax.set_xlabel(fontsize=15, rotation = 45)
    ax.set_ylabel('Percent of Occurrence', fontsize=15)
    ax.tick_params(labelsize=15)
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='darkgray')
    ax.xaxis.grid(color='darkgray')
    ax.figure.savefig("./urine_%s.pdf" %compound_name)
 
def pca(matrix, labels,features=None, k=None, body_site= None, output_loadings=None):
    print("Creating PCA Plot")
    pca = PCA(n_components=5)

    sklearn_output = pca.fit_transform(matrix)
    df = pd.DataFrame(sklearn_output)
    df["labels"] = labels
    
    
    if output_loadings is True:
        eigenvalue = (pca.explained_variance_)[0]
        eigenvalue_use = math.sqrt(eigenvalue)
        components= pca.components_
        pc1_weights = [item * eigenvalue_use for item in components[0]]
        df_temp = pd.DataFrame({"PC1_weights":pc1_weights, "Scan":features})
        df_temp.to_csv("./pc1_weights_global_individual_fecal.csv") 

    fig = plt.figure(figsize=(8.5,8.5))
    ax = fig.add_subplot()
    
    pc1 = df[0].tolist()
    pc2 = df[1].tolist()
    pc3 = df[2].tolist()
    percent_variance_explained = pca.explained_variance_ratio_
    print("Percent Variance Explained:", percent_variance_explained)
    
    plasma_cmap = matplotlib.cm.get_cmap(name="plasma")
    viridis_cmap = matplotlib.cm.get_cmap(name="viridis")

    #for numerica metadata categories we scatterplot on a gradient 
    if k == "BMI":
        labels = [float(item) for item in labels]
        plt.scatter(pc1,pc2, c=labels, cmap=viridis_cmap, s=150, edgecolors= "dimgrey")

    elif k == "perc_of_trad_diet":
        labels = [float(item) for item in labels]
        plt.scatter(pc1,pc2, c=labels, cmap=viridis_cmap, s=150, edgecolors= "dimgrey")
    
    else:
        if k == "village_of_pacient":
            cmap = ["red", "orange", "yellow", "lime", "green", "cyan", "blue", "purple", "lavender"]

        if k == "Body_site":
            cmap = ["red", "orange", "yellow", "lime", "green", "cyan", "blue", "purple", "lavender"] 
        
        if k == "Group_High_Ame_Mest":
            cmap = ["green", "blue", "orange", "purple", "magenta", "red", "yellow", "cyan", "grey"]
        
        if k in ["Cultural_individual_urbanization_groups","Global_individual_urbanization_groups", "Biological_individual_urbanization_groups", "Group_urban"]:
            cmap = URBANIZATION_COLORS

        labels = [str(item) for item in labels] #coverting type so everything plays nice
        for count,g in enumerate(np.unique(labels)):
            i = [count for count,item in enumerate(labels) if str(item) == str(g)]
            for item in i:
                if g == 'nan':
                    alpha = 0.25
                else:
                    alpha = 1
                ax.scatter(pc1[item], pc2[item],c=cmap[count], alpha=alpha, s=150, edgecolors='dimgrey')
                
    #plt.colorbar()
    ax.set_xlabel('PC1', fontsize=15)
    ax.set_ylabel('PC2', fontsize=15)
    ax.tick_params(labelsize=15)
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='darkgray')
    ax.xaxis.grid(color='darkgray')
    
    plt.show()
    
    #formatting the output files
    if body_site != None:
        ax.figure.savefig("./pca/%s_%s.pdf" %(k, body_site))    
    else:
        ax.figure.savefig("./pca/%s.pdf" %(k))    
    labels = np.unique(labels)
    num_items = len(labels) 
    colorbar_colormap(num_items)
    return(df)

def make_tree_chart():
    squarify.plot(sizes=[9,12,5,2,2], label=["High", "Mestizo", "Alto Carinagua", "Puerto Ayacucho", "Raudal de Danto"], color=[URBANIZATION_COLORS[2],URBANIZATION_COLORS[3],"red","purple", "lavender"])
    squarify.plot(sizes=[12,12], label=["Mestizo","Puerto Ayacucho"], color=[URBANIZATION_COLORS[3],"purple"])
    plt.axis('off')
    plt.show()

def make_pie_chart(matrix, labels, features):
    df = pd.read_csv("./classyfire_output.csv") 
    class_list = df["superclass"].tolist()
    df["scans"] = df["scans"].astype(str)

    unique_labels, unique_counts = np.unique(labels,return_counts=True)

    matrix = matrix.T
    reverse_df = pd.DataFrame(matrix)
    reverse_df.columns = labels
    reverse_df["scans"] = features
    reverse_df["scans"] = reverse_df["scans"].astype(str)

    combo_df = reverse_df.merge(df,left_on="scans", right_on = "scans")
    combo_df.set_index("superclass", inplace=True)

    combo_df = combo_df.T
    
    class_dict_1 = {k:0 for k in np.unique(class_list)}
    class_dict_2 = {k:0 for k in np.unique(class_list)}
    class_dict_3 = {k:0 for k in np.unique(class_list)}
    class_dict_4 = {k:0 for k in np.unique(class_list)}
    
    for index, row in combo_df.iterrows():
        for location, value in row.items():
            if value != 0:
                if index == 1:
                    class_dict_1[location] += 1
                if index == 2:
                    class_dict_2[location] += 1
                if index == 3:
                    class_dict_3[location] += 1
                if index == 4:
                    class_dict_4[location] += 1

    ploty = pd.DataFrame([class_dict_1, class_dict_2, class_dict_3, class_dict_4])
    ploty.drop("Alkaloids and derivatives", axis=1, inplace=True)

    names = ploty.columns.tolist()
    names.remove("None")

    ploty = ploty.T       
    barWidth = 0.2
    
    bars1 = [item/unique_counts[0] for item in ploty[0].tolist()]
    bars2 = [item/unique_counts[1] for item in ploty[1].tolist()]
    bars3 = [item/unique_counts[2] for item in ploty[2].tolist()]
    bars4 = [item/unique_counts[3] for item in ploty[3].tolist()]
   
    plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
    
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]

    plt.bar(r1, bars1, color=URBANIZATION_COLORS[0], width=barWidth, edgecolor='white', label='Low')
    plt.bar(r2, bars2, color=URBANIZATION_COLORS[1], width=barWidth, edgecolor='white', label='Medium')
    plt.bar(r3, bars3, color=URBANIZATION_COLORS[2], width=barWidth, edgecolor='white', label='High')
    plt.bar(r4, bars4, color=URBANIZATION_COLORS[3], width=barWidth, edgecolor='white', label='Mestizo')
      
    # Add xticks on the middle of the group bars
    plt.xlabel('Urbanization', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], names)
    
    #putting a white grid in the background 
    #plt.Axes.set_axisbelow(True)
    #ax.yaxis.grid(color='darkgray')
    #ax.xaxis.grid(color='darkgray')
         
    #rotate tick marks
    plt.xticks(rotation=15)

    #reate legend & Show graphic
    plt.legend()
    plt.show()

    #plt.savefig('urine_group.pdf')
###
def classyfire_pipeline_with_tags(quant_df):
    quant_inchi = quant_df["INCHI"].tolist()
    quant_compound_name = quant_df["Compound_Name_2"].tolist()
    quant_scans = quant_df["scans"].tolist()

    new_quant = pd.DataFrame({"Compound_Name_2":quant_compound_name, "scans":quant_scans, "INCHI":quant_inchi})

    tags = pd.read_csv("./GNPS Tag Template - MASTER - Master.csv")
    tags["Inchi_whole string"] = [str(item).replace("InChI=","") for item in tags['Inchi_whole string'].tolist()]
    tags.drop_duplicates(subset=['GNPS_annotation'], inplace=True)
    classyfire_df = new_quant.merge(tags, how="left", left_on = "Compound_Name_2", right_on='GNPS_annotation')
    og_inchi = [str(item).replace("InChI=","") for item in classyfire_df["INCHI"].tolist()] 
    tag_inchi = [str(item).replace("InChI=","") for item in classyfire_df["Inchi_whole string"].tolist()] 
    new_inchi = []

    for item1, item2 in zip(og_inchi, tag_inchi):
        if str(item1) == 'nan':
            new_inchi.append(item2)
        else:
            new_inchi.append(item1)

    classyfire_df["INCHI"]=new_inchi

    classyfire_df.to_csv("./classyfire_upload.csv")


###the main function in which we pass data/results through
def main():
    #define the parameters 
    parser = argparse.ArgumentParser(description='Do some stuff and hope it works.')
    parser.add_argument('metadata_category', help='how to color the PCA')
    parser.add_argument('output_loadings', nargs='?', default=None,help='whether or not to spit out a loadings file')
    parser.add_argument('--body_site', nargs='?', default=None)
    parser.add_argument('--fancy_stats', nargs='?', default=None)

    args = parser.parse_args()

    #parse out the parameters
    metadata_category = args.metadata_category
    body_site = args.body_site
    output_loadings = args.output_loadings
    fancy_stats = args.fancy_stats

    #filter the quant table for features that actually have an identification
    identified_feature_df = pd.read_table("./FBMN_ids.tsv")

    scans = identified_feature_df["#Scan#"].tolist()
    precollapse_compounds = identified_feature_df["Compound_Name"].tolist()
    precollapse_INCHI = identified_feature_df["INCHI"].tolist()
    precollapse_SMILES = identified_feature_df["Smiles"].tolist()

    quant_df = pd.read_csv("./GNPS_quant.csv")
    scan_df = pd.DataFrame({"scans":scans, "INCHI": precollapse_INCHI, "SMILES": precollapse_SMILES, "Compound_Name":precollapse_compounds})
    quant_df = scan_df.merge(quant_df, how="left", left_on="scans", right_on="row ID")
    
    #new feature, collapse annotations based on INCHI/SMILES/name
    quant_df = collapse_annotations(quant_df)
    quant_df.reset_index(inplace=True) 
    features = quant_df["scans"].tolist() 
    
    classyfire_pipeline_with_tags(quant_df)
    quant_df.drop(["INCHI", "Compound_Name_2", "Unnamed: 434"], inplace=True, axis=1)
    quant_df.set_index("scans",inplace=True)

    quant_df.columns = [item.replace(" Peak area","").replace(".mzXML","") for item in quant_df.columns.tolist()]
    metadata_df = pd.read_csv("./Amerindian_Metadata.csv")

    quant_cols = quant_df.columns.tolist() 
    sample_ID_quant = quant_df.columns.tolist()
    
    matrix = quant_df.to_numpy()
    matrix = matrix[:, ~np.isnan(matrix).any(axis=0)] 
    matrix = normalize_quantitative(matrix) 
    
    df = pd.DataFrame(matrix)
    df.columns=sample_ID_quant
    
    headers = (metadata_df.columns.tolist())

    #split the data according to a metadata category
    matrix, labels = split_data(df, metadata_df, metadata_category, body_site)
    labels = labels[0]
 
    #make_pie_chart(matrix, labels, features)
    #make_tree_chart()
    df = pca(matrix,labels, features, metadata_category,body_site) 
    
    if fancy_stats != None:
        #piece of code meant to retrun data for RMSSTD processing
        df.set_index('labels', inplace=True)
        unique_labels = np.unique(df.index.values.tolist())
        x_coordinates = df[0] 
        y_coordinates = df[1]
        SSW = 0 #sum of squares between groups

        RMSSTD_dict = {}
        RS_dict = {}
        removed_labels = []
        for label in unique_labels:
            cluster_x = x_coordinates[label]
            cluster_y = y_coordinates[label]

            try:
                RMSSTD_value = calculate_RMSSTD(cluster_x, cluster_y) 
                RMSSTD_dict[label] = str(RMSSTD_value)
                SSW += calculate_SSW(cluster_x, cluster_y)
            
            except:
                removed_labels.append(label)
                print("RMSSTD Calculation Failed")
        
        print(removed_labels)
        df.drop(removed_labels, axis=0, inplace=True)
        print(df)
        x_coordinates = df[0]
        y_coordinates = df[1]

        SST, x_mean, y_mean = calcaulte_SST(x_coordinates, y_coordinates)

        RS = (SST - SSW) / SST
        print(RS)
        print(RMSSTD_dict)
        
    scan = ['38560']
    compound_name = 'Fluticasone propionate'
    bar_chart_general(matrix, labels, features, scan, compound_name)

if __name__ == "__main__":
    main()
