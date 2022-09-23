import splusdata
import pandas as pd

conn = splusdata.connect("Luis", "plutarco*80") ## from splus.cloud

df = conn.query("SELECT detection.ID FROM idr3.detection_image as detection JOIN idr3.u_band as u ON detection.ID=u.ID JOIN idr3.j0378_band as J0378 ON detection.ID=J0378.ID JOIN idr3.j0395_band as J0395 ON detection.ID=J0395.ID JOIN idr3.j0410_band as J0410 ON detection.ID=J0410.ID JOIN idr3.j0430_band as J0430 ON detection.ID=J0430.ID JOIN idr3.g_band as g ON detection.ID=g.ID JOIN idr3.j0515_band as J0515 ON detection.ID=J0515.ID JOIN idr3.r_band as r ON detection.ID=r.ID JOIN idr3.j0660_band as J0660 ON detection.ID=J0660.ID JOIN idr3.i_band as i ON detection.ID=i.ID JOIN idr3.j0861_band as J0861 ON detection.ID=J0861.ID JOIN idr3.z_band as z ON detection.ID=z.ID WHERE r_PStotal < 14.0 AND e_r_PStotal <= 0.2 AND e_J0660_PStotal <= 0.2 AND e_i_PStotal <= 0.2") ## async query on splus.cloud

#print(type(df))
#df = pd.DataFrame(eval(df))
# Save the dataframe file
df.write('DR3_noFlag_3ferr_14r.csv', format='csv', overwrite=True)
#df.to_csv("DR3_noFlag_3ferr_16r.csv") 

# df = conn.query("SELECT detection.Field, detection.ID, detection.RA, detection.DEC, detection.FWHM, detection.FWHM_n, detection.A, detection.B, detection.ISOarea, detection.KRON_RADIUS, detection.nDet_PStotal, detection.PhotoFlagDet, u.u_PStotal, J0378.J0378_PStotal, J0395.J0395_PStotal, J0410.J0410_PStotal, J0430.J0430_PStotal, g.g_PStotal, J0515.J0515_PStotal, r.r_PStotal, J0660.J0660_PStotal, i.i_PStotal, J0861.J0861_PStotal, z.z_PStotal, u.e_u_PStotal, J0378.e_J0378_PStotal, J0395.e_J0395_PStotal, J0410.e_J0410_PStotal, J0430.e_J0430_PStotal, g.e_g_PStotal, J0515.e_J0515_PStotal, r.e_r_PStotal, J0660.e_J0660_PStotal, i.e_i_PStotal, J0861.e_J0861_PStotal, z.e_z_PStotal, u.u_iso, J0378.J0378_iso, J0395.J0395_iso, J0410.J0410_iso, J0430.J0430_iso, g.g_iso, J0515.J0515_iso, r.r_iso, J0660.J0660_iso, i.i_iso, J0861.J0861_iso, z.z_iso, u.e_u_iso, J0378.e_J0378_iso, J0395.e_J0395_iso, J0410.e_J0410_iso, J0430.e_J0430_iso, g.e_g_iso, J0515.e_J0515_iso, r.e_r_iso, J0660.e_J0660_iso, i.e_i_iso, J0861.e_J0861_iso, z.e_z_iso FROM idr3.detection_image as detection JOIN idr3.u_band as u ON detection.ID=u.ID JOIN idr3.j0378_band as J0378 ON detection.ID=J0378.ID JOIN idr3.j0395_band as J0395 ON detection.ID=J0395.ID JOIN idr3.j0410_band as J0410 ON detection.ID=J0410.ID JOIN idr3.j0430_band as J0430 ON detection.ID=J0430.ID JOIN idr3.g_band as g ON detection.ID=g.ID JOIN idr3.j0515_band as J0515 ON detection.ID=J0515.ID JOIN idr3.r_band as r ON detection.ID=r.ID JOIN idr3.j0660_band as J0660 ON detection.ID=J0660.ID JOIN idr3.i_band as i ON detection.ID=i.ID JOIN idr3.j0861_band as J0861 ON detection.ID=J0861.ID JOIN idr3.z_band as z ON detection.ID=z.ID WHERE r_PStotal < 14.0 AND e_r_PStotal <= 0.2 AND e_J0660_PStotal <= 0.2 AND e_i_PStotal <= 0.2") ## async query on splus.cloud