# manual entry of color
lc_color_dict = {}
for surv in ['ztf','ps','sdss']:
  lc_color_dict['g.'+surv] = 'g'
  lc_color_dict['r.'+surv] = 'r'
  lc_color_dict['i.'+surv] = 'brown'

lc_color_dict['UVW2.uvot'] = 'violet'
lc_color_dict['UVM2.uvot'] = 'magenta'
lc_color_dict['UVW1.uvot'] = 'fuchsia'

lc_color_dict['U.uvot'] = 'darkblue'
lc_color_dict['u.sdss'] = 'darkblue'

lc_color_dict['F125LP'] = 'darkviolet'
lc_color_dict['F150LP'] = 'darkviolet'
lc_color_dict['F225W'] = 'magenta'

lc_color_dict['FUV'] = 'darkviolet'
lc_color_dict['NUV'] = 'magenta'

lc_color_dict['B.uvot'] = 'lightblue'
lc_color_dict['V.uvot'] = 'orange'
lc_color_dict['c.atlas'] = 'cyan'
lc_color_dict['o.atlas'] = 'orange'

marker_dict = {key:'o' for key in lc_color_dict}
marker_dict['UVW1.uvot'] = 's'
marker_dict['UVM2.uvot'] = 's'
marker_dict['UVW2.uvot'] = 's'

marker_dict['r.ztf'] = 's'

marker_dict['F125LP'] = '*'
marker_dict['F150LP'] = 'd'
marker_dict['F225W'] = '*'
