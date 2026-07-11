'reinit'
'open tc.ctl'
'set xlopts 1 7 0.18'
'set ylopts 1 7 0.18'
'set parea 2.2 8 1.2 7'
'set mproj off'
'set grid off'
'set t 1 1'
*
* time t = 0
*
itime=1
while (itime<=1)
  'set grads off'
  'set t 'itime
  'set lat -500 0'
  'set lon 0 500'
  'set xaxis 0 1 0.2'
  'set yaxis -1 0 0.2'
  'set gxout shaded'
  'set clevs -1.0 -0.8 -0.6 -0.4 -0.2 0.2 0.4 0.6 0.8 1.0'
  'set ccols 9 14 4 11 5 0 3 10 7 12 2'
  'd p/2'
  'cbarn.gs'
  'set gxout contour'
  'set cthick 5'
  'set ccolor 1'
  'set arrscl 0.5 20'
  'set arrlab off'
  'set clevs 0'
  'set ccolor 1'
  'd lat'
  'set clevs 0'
  'set ccolor 1'
  'd lon'
  'draw xlab X (nondimensionalized)'
  'draw ylab Y (nondimensionalized)'
  'enable print out.gmf'
  'print'
  'disable print'
  '!gxps -c -i out.gmf -o out.ps'
  '!convert +antialias -rotate 90 -density 300 -geometry 1200x1200 out.ps fig3.pdf'
  itime=itime+1
*  'c'
endwhile


