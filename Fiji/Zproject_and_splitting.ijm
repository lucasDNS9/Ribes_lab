//Elie BALLOUL 2015 elie.balloul@logoss.eu

//In the middle part of the file you can change the commands
//And put new ones from the record macro command

IGet=getDirectory("Which directory to treat?");
sep=File.separator();
ISave=File.getParent(IGet) + sep + "split_channels" + File.getName(IGet) + sep; 
if(File.isDirectory(ISave)==0) File.makeDirectory(ISave);

function FindAndSave(Input, Output) {
ListFile=getFileList(Input);
for(n=0; n<ListFile.length; n++)
	{
	Target=Input + ListFile[n];
	if(File.isDirectory(Target)) {
	ITSave=Output + sep + ListFile[n] + sep; 
	if(File.isDirectory(ITSave)==0) File.makeDirectory(ITSave);
	FindAndSave(Target, ITSave);	
	}
	else {
	open(Target);
	varName=getTitle();
	OriCurrent=getImageID();


	//////////////
	///
	//You can change the next commands
title = getTitle();
run("Z Project...", "projection=[Average Intensity]");
run("Split Channels");

//Select blue canal (C=1)
selectWindow("C1-AVG_"+title);
rename("C1");
run("Red");
saveAs("Tiff", ISave + "C1-SOX2-"+title);
//waitForUser("set the color then click on OK");

//choix C=2 
selectWindow("C2-AVG_"+title);
rename("C2");
run("Green");
saveAs("Tiff", ISave + "C2-HUCD-"+title);
//waitForUser("set the color then click on OK");

//choix C=3 
selectWindow("C3-AVG_"+title);
rename("C3");
run("Blue");
saveAs("Tiff", ISave + "C3-DAPI-"+title);
//waitForUser("set the color then click on OK");

//choix C=4
//selectWindow("C4-AVG_"+title);
//rename("C4");
//run("Blue");
//saveAs("Tiff", ISave + "C4-Dapi-"+title);
//waitForUser("set the color then click on OK");

	/////// The end

	for(i=1; i<=nImages; i++) {
//	selectImage(i);
  //      ToSave=getTitle();
	//	if(ToSave==varName) {
		//	
//		}
//		else {
//	selectImage(i);
//	run("RGB Color");
//        ToSave=getTitle();
 //       save(Output + ToSave);
//		}
	}
	close("*");
	}
	}
}
FindAndSave(IGet, ISave);