from pathlib import Path
import numpy as np,json
from skimage.measure import regionprops
P=Path(__file__).resolve().parent;Q=P/'qa';Q.mkdir(exist_ok=True)
def main():
    yy,xx=np.mgrid[-90:91,-90:91];out=[]
    for row_angle in [0.,.3,.6,-.5,1.2,np.pi/2]:
        along=yy*np.cos(row_angle)+xx*np.sin(row_angle);across=-yy*np.sin(row_angle)+xx*np.cos(row_angle);mask=(along/60)**2+(across/17)**2<=1;theta=regionprops(mask.astype(int))[0].orientation;image_angle=np.pi/2-theta;truth=np.pi/2-row_angle;agreement=float(np.cos(2*(image_angle-truth)));assert agreement>.99999;out.append({'true_row_angle':row_angle,'regionprops_row_angle':float(theta),'converted_image_x_angle':float(image_angle),'axis_agreement':agreement})
    result={'all_pass':True,'synthetic_ellipses':out,'minimum_axis_agreement':min(r['axis_agreement'] for r in out),'scope':'Frame conversion for new positional metric and raw-image overlays; original pairwise nuclear alignment is invariant and unchanged.'};(Q/'angle_convention.json').write_text(json.dumps(result,indent=2));print('Angle-frame synthetic checks passed',result['minimum_axis_agreement'])
if __name__=='__main__':main()
