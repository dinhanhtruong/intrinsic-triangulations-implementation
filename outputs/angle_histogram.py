import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import argparse

parser = argparse.ArgumentParser(description="Converts the raw txt corner angle files into histogram images. The second image will contain the overlaid histograms")
parser.add_argument("img_name", type=str, help="name of the output mesh file without the extension (e.g. 'bread')")
args = parser.parse_args()

for i in range(2):
    mesh_type = "extrinsic" if i == 0 else "intrinsic"
    # load angles 
    angles = np.loadtxt(f"{args.img_name}_{mesh_type}_angles.txt")
    n = len(angles)
    # convert to degrees
    angles *= 360/(2*np.pi)

    plt.hist(angles, bins=18, range=[0, 180], weights=np.ones(n)/n, alpha=0.7) # normalize frequencies
    plt.xticks(np.arange(0, 181, 20))
    # format the y-axis to display percentage of total num angles
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    plt.title("Corner Angles")
    plt.xlabel("Angle (degrees)")
    plt.ylabel("Relative Frequency")

    save_path = f'{args.img_name}_{mesh_type}_angles.png'
    plt.savefig(save_path, dpi=300)
    print(f"Saved histogram at {save_path}")
    # plt.show()