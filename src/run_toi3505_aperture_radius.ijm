argument = getArgument();
settings = split(argument, ";");
if (settings.length != 5)
    exit("Use: image_directory;source_radius;outer_sky_radius;output_table;status_file");

imageDirectory = settings[0];
sourceRadius = parseInt(settings[1]);
outerSkyRadius = parseInt(settings[2]);
outputTable = settings[3];
statusFile = settings[4];

File.saveString(
    "state=starting\n" +
    "source_radius=" + sourceRadius + "\n" +
    "sky_ring=70-" + outerSkyRadius + "\n",
    statusFile
);

run("Image Sequence...", "open=[" + imageDirectory + "] sort use");
if (!is("Virtual Stack"))
    exit("The images did not open as a virtual stack.");

call("ij.Prefs.set", "aperture.radius", "" + sourceRadius);
call("ij.Prefs.set", "aperture.rback1", "70");
call("ij.Prefs.set", "aperture.rback2", "" + outerSkyRadius);
call("ij.Prefs.set", "aperture.ccdgain", "1.0");
call("ij.Prefs.set", "aperture.ccdnoise", "1.414214");
call("ij.Prefs.set", "aperture.ccddark", "0.012283");
call("ij.Prefs.set", "aperture.darkkeyword", "");
call("ij.Prefs.set", "aperture.reposition", "true");
call("ij.Prefs.set", "aperture.removebackstars", "true");
call("ij.Prefs.set", "aperture.backplane", "false");
call("ij.Prefs.set", "aperture.showremovedpixels", "false");
call("ij.Prefs.set", "aperture.fitskeywords", "BJD_TDB,AIRMASS");
call("ij.Prefs.set", "aperture.showtimes", "true");
call("ij.Prefs.set", "aperture.showfits", "true");
call("ij.Prefs.set", "aperture.showposition", "true");
call("ij.Prefs.set", "aperture.showpositionfits", "true");
call("ij.Prefs.set", "aperture.showphotometry", "true");
call("ij.Prefs.set", "aperture.showback", "true");
call("ij.Prefs.set", "aperture.showwidths", "true");
call("ij.Prefs.set", "aperture.showmeanwidth", "true");
call("ij.Prefs.set", "aperture.showerrors", "true");
call("ij.Prefs.set", "aperture.showsnr", "true");

call("ij.Prefs.set", "multiaperture.finished", "false");
call("ij.Prefs.set", "multiaperture.canceled", "false");
call("ij.Prefs.set", "multiaperture.automode", "true");
call("ij.Prefs.set", "multiaperture.previous", "false");
call("ij.Prefs.set", "multiaperture.singlestep", "false");
call("ij.Prefs.set", "multiaperture.allowsinglestepapchanges", "false");
call("ij.Prefs.set", "multiaperture.usevarsizeap", "false");
call("ij.Prefs.set", "multiaperture.apradius", "FIXED");
call("ij.Prefs.set", "multiaperture.usema", "true");
call("ij.Prefs.set", "multiaperture.usealign", "false");
call("ij.Prefs.set", "multiaperture.usewcs", "false");
call("ij.Prefs.set", "multiaperture.haltOnError", "true");
call("ij.Prefs.set", "multiaperture.suggestCompStars", "false");
call("ij.Prefs.set", "multiaperture.showhelp", "false");
call("ij.Prefs.set", "multiaperture.getMags", "false");
call("ij.Prefs.set", "multiaperture.updatePlot", "false");
call("ij.Prefs.set", "multiaperture.naperturesmax", "11");

call("ij.Prefs.set", "multiaperture.xapertures", "1850.080000,2719.759251,1698.820436,658.556838,3131.451150,2078.705977,3331.827655,636.135202,2003.188107,2301.668841,2514.389245");
call("ij.Prefs.set", "multiaperture.yapertures", "2353.680000,929.371245,3526.551310,790.088175,2253.832657,843.098419,3590.876370,1726.370761,1316.363418,3044.178143,2131.728407");
call("ij.Prefs.set", "multiaperture.isrefstar", "false,true,true,true,true,true,true,true,true,true,true");
call("ij.Prefs.set", "multiaperture.isalignstar", "true,true,true,true,true,true,true,true,true,true,true");
call("ij.Prefs.set", "multiaperture.centroidstar", "true,true,true,true,true,true,true,true,true,true,true");
call("ij.Prefs.set", "multiaperture.raapertures", "");
call("ij.Prefs.set", "multiaperture.decapertures", "");
call("ij.Prefs.set", "multiaperture.absmagapertures", "");

File.append("state=measuring\n", statusFile);
run("MultiAperture ");
elapsedSeconds = 0;
wait(1000);
selectWindow("Measurements");
while (Table.size < 281) {
    wait(1000);
    elapsedSeconds++;
    if (elapsedSeconds > 900)
        exit("Multi-Aperture did not finish within 15 minutes.");
}
File.append("state=measurement_complete\nelapsed_seconds=" + elapsedSeconds + "\n", statusFile);
Table.save(outputTable);
File.append("state=table_saved\n", statusFile);
