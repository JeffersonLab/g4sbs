/**
 * @file G4SBSGrinch.cc
 * @brief Source codes for detector construction
 * @author Huan Yao <hyao@jlab.org>
 * @version 1.0
 * @date 2012-05-04
 */
#include "G4SBSGrinch.hh"

#define CHAR_LEN 255

enum E_VERTICES {
	k_BOTTOM_LEFT,
	k_TOP_LEFT,
	k_TOP_RIGHT,
	k_BOTTOM_RIGHT
};

/* --------------------------------------------------------------------------*/
/**
 * @brief Calculate all vertices on two trapezoid planes used in BB_EM and BB_M
 *
 * @param aTop Top Box full size
 * @param aBottom Bottom Box full size
 * @param aEntrance Entrance Box full size
 * @param aCorner_Bend_Angle Corner Bending angle
 * @param aBB_Bottom_Length BB Bottom Length defined in BGCUserManager.hh
 *
 * @return vector<G4TwoVector> output
 */
inline std::vector<G4TwoVector> Two_Trapezoid_Planes(const G4ThreeVector& aTop,const G4ThreeVector& aBottom,const G4ThreeVector& aEntrance,
		const G4double& aCorner_Bend_Angle,const G4double& aBB_Bottom_Length) {
    std::vector<G4TwoVector> output; //0=4, 1=5, 2=6, 3=7
	G4TwoVector vertices;
	G4int i;
	for ( i = 0; i < 2; ++i ) { //2=two planes
		//0/4
		vertices.set(-aBB_Bottom_Length*0.5,-aEntrance.y()*0.5);
		output.push_back(vertices);
		//1/5
		vertices.set(-aBB_Bottom_Length*0.5,aEntrance.y()*0.5);
		output.push_back(vertices);
		//2/6
		vertices.set(aTop.x()*cos(aCorner_Bend_Angle)-aBB_Bottom_Length*0.5,aTop.x()*sin(aCorner_Bend_Angle)+aEntrance.y()*0.5);
		output.push_back(vertices);
		//3/7
		vertices.set(aBottom.x()-aBB_Bottom_Length*0.5,-aEntrance.y()*0.5);
		output.push_back(vertices);
	}

	return output;
}

inline G4VSolid* Construct_GC_Tank_Box(G4String aName, G4ThreeVector aInner_Full_Size, G4double aThickness, G4ThreeVector aEntrance_Window_Full_Size, G4ThreeVector aExit_Window_Full_Size)
{
	G4Box*                  box;
	G4Box*                  cutbox;
	G4Box*                  cutentrance;
	G4Box*                  cutexit;
	G4SubtractionSolid*     sub_solid0;
	G4SubtractionSolid*     sub_solid1;
	G4SubtractionSolid*     sub_solid2;

	G4ThreeVector trans;
	G4RotationMatrix tempRM = tempRM.IDENTITY;

	box=new G4Box("box",aInner_Full_Size.x()*0.5+aThickness,aInner_Full_Size.y()*0.5+aThickness,aInner_Full_Size.z()*0.5+aThickness);
	cutbox=new G4Box("cutebox",aInner_Full_Size.x()*0.5,aInner_Full_Size.y()*0.5,aInner_Full_Size.z()*0.5);
	cutentrance=new G4Box("cutentrance",(aEntrance_Window_Full_Size.x()+aThickness)*0.5,aEntrance_Window_Full_Size.y()*0.5,aEntrance_Window_Full_Size.z()*0.5);
	cutexit=new G4Box("cutexit",(aExit_Window_Full_Size.x()+aThickness)*0.5,aExit_Window_Full_Size.y()*0.5,aExit_Window_Full_Size.z()*0.5);

	sub_solid0 = new G4SubtractionSolid("sub_solid0", box, cutbox);

	trans.set(-aInner_Full_Size.x()*0.5-aThickness*0.5,0,0);
	sub_solid1 = new G4SubtractionSolid("sub_solid1", sub_solid0, cutentrance, G4Transform3D(tempRM,trans));

	trans.set(aInner_Full_Size.x()*0.5+aThickness*0.5,0,0);
	sub_solid2 = new G4SubtractionSolid(aName, sub_solid1, cutexit, G4Transform3D(tempRM,trans));

	return sub_solid2;
}

/* --------------------------------------------------------------------------*/
/**
 * @brief Generate a mirror at (diameter*0.5,0,0) with height=Vert_size=z-axis, Width=hori_size=y-axis, -x-axis is normal to surface of mirror and x=radius_mirror=diameter*0.5
 * So need to rotate mirror about x-axis -90 deg to make mirror at right position in real life
 * Do not use G4Tubs directly just because I want the edge to be parallel to vertical/horizontal axis
 * @param sName Mirror Name
 * @param diameter diameter of mirror
 * @param vert_size vertical full size
 * @param hori_size horizontal full size
 * @param thickness thickness of mirror
 *
 * @return G4VSolid
 */
inline G4VSolid* Construct_CylinderMirror(G4String sName, G4double diameter, G4double vert_size, G4double hori_size, G4double thickness)
{
	G4Box*                  cutbox;
	G4Tubs*                 tube;
	G4SubtractionSolid*     sub_solid0;
	G4SubtractionSolid*     sub_solid1;
	G4SubtractionSolid*     sub_solid2;
	G4SubtractionSolid*     mirror_solid;

	G4ThreeVector trans;
	G4RotationMatrix tempRM = tempRM.IDENTITY;

	tube = new G4Tubs("tube", diameter*0.5, diameter*0.5+thickness, vert_size*0.5, -90*deg, 180*deg);
	cutbox = new G4Box("cutbox", diameter, diameter, diameter);

	trans.set(0, diameter+hori_size*0.5,0);
	sub_solid0 = new G4SubtractionSolid("sub_solid0", tube, cutbox, G4Transform3D(tempRM,trans));

	trans.set(0, -(diameter+hori_size*0.5),0);
	sub_solid1 = new G4SubtractionSolid("sub_solid1", sub_solid0, cutbox, G4Transform3D(tempRM,trans));

	trans.set(0, 0,diameter+vert_size*0.5);
	sub_solid2 = new G4SubtractionSolid("sub_solid2", sub_solid1, cutbox, G4Transform3D(tempRM,trans));

	trans.set(0, 0,-(diameter+vert_size*0.5));
	mirror_solid = new G4SubtractionSolid(sName, sub_solid2, cutbox, G4Transform3D(tempRM,trans));

	return mirror_solid;
}

inline G4VSolid* Construct_SphereMirror(G4String sName, G4double diameter, G4double vert_size, G4double hori_size, G4double thickness, G4double startDx)
{
	G4Box*                  cutbox;
	G4Sphere*               sphere;
	G4SubtractionSolid*     sub_solid0;
	G4SubtractionSolid*     sub_solid1;
	G4SubtractionSolid*     sub_solid2;
	G4SubtractionSolid*     mirror_solid;

	G4ThreeVector trans;
	G4RotationMatrix tempRM = tempRM.IDENTITY;

	sphere = new G4Sphere("sphere", diameter*0.5, diameter*0.5+thickness, 0, 180*deg, 0*deg, 180*deg);
	cutbox = new G4Box("cutbox", diameter, diameter, diameter);

	if(vert_size > 0)
	{
		trans.set(0, 0, -diameter);
		sub_solid0 = new G4SubtractionSolid("sub_solid0", sphere, cutbox, G4Transform3D(tempRM,trans));

		trans.set(0, 0, (vert_size+diameter));
		sub_solid1 = new G4SubtractionSolid("sub_solid1", sub_solid0, cutbox, G4Transform3D(tempRM,trans));

		trans.set(diameter+startDx, 0, 0);
		sub_solid2 = new G4SubtractionSolid("sub_solid2", sub_solid1, cutbox, G4Transform3D(tempRM,trans));

		trans.set(-diameter+startDx-hori_size, 0, 0);
		mirror_solid = new G4SubtractionSolid(sName, sub_solid2, cutbox, G4Transform3D(tempRM,trans));
	}

	if(vert_size < 0)
	{
		trans.set(0, 0, diameter);
		sub_solid0 = new G4SubtractionSolid("sub_solid0", sphere, cutbox, G4Transform3D(tempRM,trans));

		trans.set(0, 0, (vert_size-diameter));
		sub_solid1 = new G4SubtractionSolid("sub_solid1", sub_solid0, cutbox, G4Transform3D(tempRM,trans));

		trans.set(diameter+startDx, 0, 0);
		sub_solid2 = new G4SubtractionSolid("sub_solid2", sub_solid1, cutbox, G4Transform3D(tempRM,trans));

		trans.set(-diameter+startDx-hori_size, 0, 0);
		mirror_solid = new G4SubtractionSolid(sName, sub_solid2, cutbox, G4Transform3D(tempRM,trans));
	}

	return mirror_solid;

}

/* --------------------------------------------------------------------------*/
/**
 * @brief Construct square-opening cones
 *
 * @param aName Name of solid
 * @param aRmin1 inner radius of smaller opening at -aDz
 * @param aRmax1 outer radius of smaller opening at -aDz
 * @param aRmin2 inner radius of large opening at aDz
 * @param aRmax2 outer radius of large opening at aDz
 * @param aDz half length of cone
 * @param aSPhi start phi angle of cone
 * @param aDPhi phi angle range
 * @param aVSize vertical size between each row of pmt
 * @param aHSize hori_size between two adjacent pmts
 *
 * @return G4VSolid
 */
inline G4VSolid* Construct_Square_Opening_Cone(G4String aName, G4double aRmin1, G4double aRmax1, G4double aRmin2, G4double aRmax2, G4double aDz, G4double aSPhi, G4double aDPhi, G4double aVSize, G4double aHSize)
{
	G4Box*                  box;
	G4Cons*                 cone;
	G4SubtractionSolid*     sub_solid0;

	G4ThreeVector trans;
	G4RotationMatrix tempRM = tempRM.IDENTITY;

	//cone=new G4Cons("cone",aRmin1,aRmax1,aRmin2,aRmax2,aDz,aSPhi,aDPhi);
	//box = new G4Box("box", aHSize*0.5, aVSize*0.5, aDz+1);

	//trans.set(aHSize,0,0);
	//sub_solid0 = new G4SubtractionSolid("sub_solid0", cone, box, G4Transform3D(tempRM,trans));

	//trans.set(-aHSize,0,0);
	//sub_solid1 = new G4SubtractionSolid("sub_solid1", sub_solid0, box, G4Transform3D(tempRM,trans));

	//trans.set(0,aVSize,0);
	//sub_solid2 = new G4SubtractionSolid("sub_solid2", sub_solid1, box, G4Transform3D(tempRM,trans));

	//trans.set(0,-aVSize,0);
	//cone_solid = new G4SubtractionSolid(aName, sub_solid2, box, G4Transform3D(tempRM,trans));
	//
	//+1 to cut box thoroughly,otherwise the outer surface will still exist
	G4double extended_Rmin1=((2*aDz+1)*aRmin1-aRmin2)/(2*aDz);
	G4double extended_Rmin2=((2*aDz+1)*aRmin2-aRmin1)/(2*aDz);
	G4double extended_Rmax1=((2*aDz+1)*aRmax1-aRmax2)/(2*aDz);
	G4double extended_Rmax2=((2*aDz+1)*aRmax2-aRmax1)/(2*aDz);

	cone=new G4Cons("cone",extended_Rmin1,extended_Rmax1,extended_Rmin2,extended_Rmax2,aDz+1,aSPhi,aDPhi);
	box = new G4Box("box", aHSize*0.5, aVSize*0.5, aDz);

	trans.set(0,0,0);
	sub_solid0 = new G4SubtractionSolid("sub_solid0", box, cone, G4Transform3D(tempRM,trans));

	return sub_solid0;
}


/* --------------------------------------------------------------------------*/
/**
 * @brief Construct Function
 */
G4SBSGrinch::G4SBSGrinch(G4SBSDetectorConstruction *dc):G4SBSComponent(dc) {
}

G4SBSGrinch::~G4SBSGrinch() {
}

/* --------------------------------------------------------------------------*/
/**
 * @brief Main construct method
 *  Center of world is the center of bigbite magnet
 * @return Always world physics
 */
void  G4SBSGrinch::BuildComponent(G4LogicalVolume *bblog) {

	G4ThreeVector Translation(0,0,0);
	G4RotationMatrix Rotation;
	G4int i,j;



	G4RotationMatrix rm, PMTrm;
	G4Transform3D movetocenter,offset,combine;


	char tmpname[CHAR_LEN];


	/* --------------------------------------------------------------------------*/
	//GC_Tank
	G4String GC_Tank_Name("GC_Tank");
	G4String GC_Tank_Material=("C4F8O");
	G4ThreeVector GC_Tank_Inner_FullSize(92*cm, 200*cm, 150.229*cm);
	G4double GC_Tank_Thickness= 1.27*cm;
	G4VSolid* GC_Tank_Solid=ConstructSimple(GC_Tank_Name,G4String("G4Box"),GC_Tank_Inner_FullSize);
	GC_Tank_log=new G4LogicalVolume(GC_Tank_Solid,GetMaterial(GC_Tank_Material),GC_Tank_Name+"_log");

	G4VisAttributes* GC_Tank_log_VisAtt = new G4VisAttributes();
	GC_Tank_log_VisAtt->SetColor(GetColor(G4String("Cyan")));
	GC_Tank_log_VisAtt->SetVisibility(true);
	GC_Tank_log_VisAtt->SetForceWireframe(true);
	GC_Tank_log->SetVisAttributes(GC_Tank_log_VisAtt);


	/* --------------------------------------------------------------------------*/
	//GC_Tank_Entrance_Window
	G4String GC_Tank_Entrance_Window_Name("GC_Tank_Entrance_Window");
	G4String GC_Tank_Entrance_Window_Material("Mylar");
	G4ThreeVector GC_Tank_Entrance_Window_FullSize(0.01*cm, 170*cm, 45*cm);
	G4VSolid* GC_Tank_Entrance_Window_Solid=ConstructSimple(GC_Tank_Entrance_Window_Name,G4String("G4Box"),GC_Tank_Entrance_Window_FullSize);
	G4LogicalVolume* GC_Tank_Entrance_Window_log=new G4LogicalVolume(GC_Tank_Entrance_Window_Solid,GetMaterial(GC_Tank_Entrance_Window_Material),GC_Tank_Entrance_Window_Name+"_log");

	G4VisAttributes* GC_Tank_Entrance_Window_log_VisAtt = new G4VisAttributes();
	GC_Tank_Entrance_Window_log_VisAtt->SetColor(GetColor(G4String("Black")));
	GC_Tank_Entrance_Window_log_VisAtt->SetVisibility(true);
	GC_Tank_Entrance_Window_log_VisAtt->SetForceWireframe(false);
	GC_Tank_Entrance_Window_log->SetVisAttributes(GC_Tank_Entrance_Window_log_VisAtt);



	/* --------------------------------------------------------------------------*/
	//GC_Tank_Exit_Window
	G4String GC_Tank_Exit_Window_Name("GC_Tank_Exit_Window");
	G4String GC_Tank_Exit_Window_Material("Mylar");
	G4ThreeVector GC_Tank_Exit_Window_FullSize(0.01*cm, 200*cm, 55*cm);
	G4VSolid* GC_Tank_Exit_Window_Solid=ConstructSimple(GC_Tank_Exit_Window_Name,G4String("G4Box"),GC_Tank_Exit_Window_FullSize);
	G4LogicalVolume* GC_Tank_Exit_Window_log=new G4LogicalVolume(GC_Tank_Exit_Window_Solid,GetMaterial(GC_Tank_Exit_Window_Material),GC_Tank_Exit_Window_Name+"_log");

	G4VisAttributes* GC_Tank_Exit_Window_log_VisAtt = new G4VisAttributes();
	GC_Tank_Exit_Window_log_VisAtt->SetColor(GetColor(G4String("Black")));
	GC_Tank_Exit_Window_log_VisAtt->SetVisibility(true);
	GC_Tank_Exit_Window_log_VisAtt->SetForceWireframe(false);
	GC_Tank_Exit_Window_log->SetVisAttributes(GC_Tank_Exit_Window_log_VisAtt);


	//GC_Tank_Box
	G4VSolid* GC_Tank_Box_Solid=Construct_GC_Tank_Box(GC_Tank_Name+"_Box", GC_Tank_Inner_FullSize, GC_Tank_Thickness, GC_Tank_Entrance_Window_FullSize, GC_Tank_Exit_Window_FullSize);
	G4LogicalVolume* GC_Tank_Box_log=new G4LogicalVolume(GC_Tank_Box_Solid,GetMaterial("Stainless_Steel"),GC_Tank_Name+"_Box_log");

	G4VisAttributes* GC_Tank_Box_log_VisAtt = new G4VisAttributes();
	GC_Tank_Box_log_VisAtt->SetColor(GetColor(G4String("Cyan")));
	GC_Tank_Box_log_VisAtt->SetVisibility(true);
	GC_Tank_Box_log_VisAtt->SetForceWireframe(true);
	GC_Tank_Box_log->SetVisAttributes(GC_Tank_Box_log_VisAtt);



	G4AssemblyVolume* assemblyWindow = new G4AssemblyVolume();
	rm = rm.IDENTITY;


	Translation.set(0,0,0);
	assemblyWindow->AddPlacedVolume(GC_Tank_Box_log, Translation, &rm);
	Translation.set(-GC_Tank_Inner_FullSize.x()/2,0,0);
	assemblyWindow->AddPlacedVolume(GC_Tank_Entrance_Window_log, Translation, &rm);
	Translation.set(GC_Tank_Inner_FullSize.x()/2,0,0);
	assemblyWindow->AddPlacedVolume(GC_Tank_Exit_Window_log, Translation, &rm);

	Translation.set(0,0,0);
	assemblyWindow->MakeImprint(GC_Tank_log, Translation, &rm);
//	usermanager->Inc_PMT_ID_Offset(1);//count GC_Tank_Box
//	usermanager->Inc_PMT_ID_Offset(2);//count GC_Tank_Entrance_Window,GC_Tank_Exit_Window
//	usermanager->Inc_PMT_ID_Offset(1);//count GC_Tank

	/* --------------------------------------------------------------------------*/
	//GC_Mirror
	const G4int num_of_mirrors=4;
	G4VSolid** GC_Mirror=new G4VSolid*[num_of_mirrors];
	G4LogicalVolume** GC_Mirror_log=new G4LogicalVolume*[num_of_mirrors];
	G4VisAttributes** GC_Mirror_log_VisAtt=new G4VisAttributes*[num_of_mirrors];
	G4AssemblyVolume* assemblyMir = new G4AssemblyVolume();
	G4double* GC_Mirror_Radius=new G4double[num_of_mirrors];
	G4double* GC_Mirror_Ver_Size=new G4double[num_of_mirrors];
	G4double* GC_Mirror_Rot_Y_Angle=new G4double[num_of_mirrors];
	G4double* GC_Mirror_Rot_Axis_Angle=new G4double[num_of_mirrors];

	G4int mirror_index;
	G4ThreeVector newAxis;
	G4String mirror_name;

        char tempstr[255];

        double mirror_radius[num_of_mirrors] = { 130*cm, 130*cm, 130*cm, 130*cm };
        double mirror_ver_size[num_of_mirrors] = { 40*cm, 60*cm, 60*cm, 40*cm };
        double mirror_hor_size[num_of_mirrors] = { 70*cm, 70*cm, 70*cm, 70*cm };
        double mirror_thickness[num_of_mirrors] = { 0.635*cm, 0.635*cm, 0.635*cm, 0.635*cm };
        double mirror_rot_y_angle[num_of_mirrors] = { 27.5*deg, 27.5*deg, 27.5*deg, 27.5*deg };
        double mirror_rot_axis_angle[num_of_mirrors] = { 10*deg, 0*deg, 0*deg, -10*deg };
        G4String mirror_color[num_of_mirrors] = {G4String("Yellow"), G4String("Blue"), G4String("Red"), G4String("Blue")};
        G4ThreeVector mirror_new_offset[num_of_mirrors] = {G4ThreeVector(16.9194*cm, 80*cm, 1.60364*cm), G4ThreeVector(), G4ThreeVector(), G4ThreeVector(16.9194*cm, -80*cm, 1.60364*cm)  };

	for ( i = 0; i < num_of_mirrors; ++i ) {
		mirror_index=i+1;
		GC_Mirror_Radius[i]= mirror_radius[mirror_index-1];
		GC_Mirror_Ver_Size[i]= mirror_ver_size[mirror_index-1];
		GC_Mirror_Rot_Y_Angle[i]=mirror_rot_y_angle[mirror_index-1];
		GC_Mirror_Rot_Axis_Angle[i]=mirror_rot_axis_angle[mirror_index-1];
                sprintf(tempstr, "GC_Mirror_%d", mirror_index);
		mirror_name=G4String(tempstr);
		GC_Mirror[i]=Construct_CylinderMirror(mirror_name,2*GC_Mirror_Radius[i],GC_Mirror_Ver_Size[i], mirror_hor_size[mirror_index-1],mirror_thickness[mirror_index-1]);
		GC_Mirror_log[i]=new G4LogicalVolume(GC_Mirror[i], GetMaterial(G4String("Glass")), mirror_name+"_log", 0, 0, 0);
		GC_Mirror_log_VisAtt[i]=new G4VisAttributes();
		GC_Mirror_log_VisAtt[i]->SetColor(GetColor(mirror_color[mirror_index-1]));
		GC_Mirror_log_VisAtt[i]->SetVisibility(true);
		GC_Mirror_log[i]->SetVisAttributes(GC_Mirror_log_VisAtt[i]);

		rm=rm.IDENTITY;
		rm.rotateX(-90*deg);
		Translation.set(-GC_Mirror_Radius[i],0,0);
		movetocenter=G4Transform3D(rm,Translation);//After movetocenter, mirror center location is at center of tank
		rm=rm.IDENTITY;
		rm.rotateY(GC_Mirror_Rot_Y_Angle[i]);
		newAxis.set(0,0,1);//First Z-axis
		newAxis.rotateY(GC_Mirror_Rot_Y_Angle[i]);
		rm.rotate(GC_Mirror_Rot_Axis_Angle[i],newAxis);
		Translation=mirror_new_offset[mirror_index-1];
		offset=G4Transform3D(rm,Translation);//After offset, mirror center is at right position
		combine=offset*movetocenter;
		assemblyMir->AddPlacedVolume(GC_Mirror_log[i], combine);
	}
	rm = rm.IDENTITY;
	Translation.set(0,0,0);
	assemblyMir->MakeImprint(GC_Tank_log, Translation, &rm);
//	usermanager->Inc_PMT_ID_Offset(num_of_mirrors);//Count num_of_mirrors
	/* --------------------------------------------------------------------------*/



	/* --------------------------------------------------------------------------*/
	//GC PMT Box and PMT array
	G4ThreeVector GC_PMT_Box_FullSize(18.75*cm, 192*cm, 33.9*cm);
	G4double GC_PMT_Box_Thickness=0.635*cm;
	G4String GC_PMT_Box_Name("GC_PMT_Box");
	G4Box* GC_PMT_Outer_Box = new G4Box("GC_PMT_Outer_Box", (GC_PMT_Box_FullSize.x()+GC_PMT_Box_Thickness)*0.5, (GC_PMT_Box_FullSize.y()+2*GC_PMT_Box_Thickness)*0.5, (GC_PMT_Box_FullSize.z()+2*GC_PMT_Box_Thickness)*0.5);
	G4Box* GC_PMT_Inner_Box = new G4Box("GC_PMT_Inner_Box", GC_PMT_Box_FullSize.x()*0.5, GC_PMT_Box_FullSize.y()*0.5, GC_PMT_Box_FullSize.z()*0.5);
	rm = rm.IDENTITY;
	Translation.set(GC_PMT_Box_Thickness*0.5,0,0);
	G4SubtractionSolid* GC_PMT_Box = new G4SubtractionSolid(GC_PMT_Box_Name.data(), GC_PMT_Outer_Box, GC_PMT_Inner_Box, G4Transform3D(rm,Translation));
	//G4Box* GC_PMT_Box = new G4Box("GC_PMT_Box", GC_PMT_Box_FullSize.x()*0.5, GC_PMT_Box_FullSize.y()*0.5, GC_PMT_Box_FullSize.z()*0.5);
	G4LogicalVolume* GC_PMT_Box_log = new G4LogicalVolume(GC_PMT_Box, GetMaterial(G4String("Stainless_Steel")), GC_PMT_Box_Name+"_log", 0, 0, 0);
	G4VisAttributes* GC_PMT_Box_log_VisAtt = new G4VisAttributes();
	GC_PMT_Box_log_VisAtt->SetColor(GetColor(G4String("Grey")));
	GC_PMT_Box_log_VisAtt->SetVisibility(true);
	GC_PMT_Box_log_VisAtt->SetForceWireframe(true);
	GC_PMT_Box_log->SetVisAttributes(GC_PMT_Box_log_VisAtt);

	//GC PMT Jail
	G4ThreeVector GC_PMT_Box_Jail_FullSize(10*cm, 0.1*cm, 33.9*cm);
	G4double GC_PMT_Box_Jail_Space=-5.0*cm;
	G4String GC_PMT_Box_Jail_Name("GC_PMT_Box_Jail");
	G4int Num_Of_GC_PMT_Box_Jail_Bars=61;
	G4Box* GC_PMT_Box_Jail = new G4Box(GC_PMT_Box_Jail_Name.data(), GC_PMT_Box_Jail_FullSize.x()*0.5, GC_PMT_Box_Jail_FullSize.y()*0.5, GC_PMT_Box_Jail_FullSize.z()*0.5);
	G4LogicalVolume* GC_PMT_Box_Jail_log = new G4LogicalVolume(GC_PMT_Box_Jail, GetMaterial(G4String("mu-metal")), GC_PMT_Box_Jail_Name+"_log", 0, 0, 0);
	G4VisAttributes* GC_PMT_Box_Jail_log_VisAtt = new G4VisAttributes();
	GC_PMT_Box_Jail_log_VisAtt->SetColor(GetColor(G4String("Black")));
	GC_PMT_Box_Jail_log_VisAtt->SetVisibility(true);
	GC_PMT_Box_Jail_log_VisAtt->SetForceWireframe(false);
	GC_PMT_Box_Jail_log->SetVisAttributes(GC_PMT_Box_Jail_log_VisAtt);

	//GC PMT
	G4double GC_PMT_Radius = 1.25*cm;
	G4double GC_PMT_Cover_Radius = 1.5*cm;
	G4double GC_PMT_Length = 13.75*cm;
	G4String GC_PMT_Name("GC_PMT");
	G4String GC_PMT_Glass_Name=GC_PMT_Name+"_Glass";
	G4double GC_PMT_Glass_Thickness=0.3*cm;

	G4Tubs* GC_PMT = new G4Tubs(GC_PMT_Name.data(), 0, GC_PMT_Radius, (GC_PMT_Length-GC_PMT_Glass_Thickness)*0.5, 0, 360*deg);
	G4LogicalVolume* GC_PMT_log = new G4LogicalVolume(GC_PMT, GetMaterial(G4String("Al")), GC_PMT_Name+"_log", 0, 0, 0);
	G4VisAttributes* GC_PMT_log_VisAtt = new G4VisAttributes();
	GC_PMT_log_VisAtt->SetColor(GetColor(G4String("Blue")));
	GC_PMT_log_VisAtt->SetVisibility(true);
	GC_PMT_log->SetVisAttributes(GC_PMT_log_VisAtt);

	G4Tubs* GC_PMT_Cover = new G4Tubs("GC_PMT_Cover", GC_PMT_Radius, GC_PMT_Cover_Radius, GC_PMT_Length*0.5, 0, 360*deg);
	G4LogicalVolume* GC_PMT_Cover_log = new G4LogicalVolume(GC_PMT_Cover, GetMaterial("Stainless_Steel"), "GC_PMT_Cover_log", 0, 0, 0);
	G4VisAttributes* GC_PMT_Cover_log_VisAtt = new G4VisAttributes();
	GC_PMT_Cover_log_VisAtt->SetColor(GetColor(G4String("Grey")));
	GC_PMT_Cover_log_VisAtt->SetVisibility(true);
	GC_PMT_Cover_log->SetVisAttributes(GC_PMT_Cover_log_VisAtt);

	G4Tubs* GC_PMT_Glass;
	G4LogicalVolume* GC_PMT_Glass_log;
	G4VisAttributes* GC_PMT_Glass_log_VisAtt;
	if ( fabs(GC_PMT_Glass_Thickness)>1e-5 ) {
		GC_PMT_Glass = new G4Tubs(GC_PMT_Glass_Name.data(), 0, GC_PMT_Radius, GC_PMT_Glass_Thickness*0.5, 0, 360*deg);
		GC_PMT_Glass_log = new G4LogicalVolume(GC_PMT_Glass, GetMaterial(G4String("Quartz")), GC_PMT_Glass_Name+"_log", 0, 0, 0);
		GC_PMT_Glass_log_VisAtt = new G4VisAttributes();
		GC_PMT_Glass_log_VisAtt->SetColor(GetColor(G4String("Cyan")));
		GC_PMT_Glass_log_VisAtt->SetVisibility(true);
		GC_PMT_Glass_log->SetVisAttributes(GC_PMT_Glass_log_VisAtt);
	}


	//GC_PMT_Cone
	G4String GC_PMT_Cone_Name("GC_PMT_Cone");
	G4String GC_PMT_Cone_Shape("G4Cons");
	G4double GC_PMT_Cone_Max_Inner_Radius=2.15555*cm;
	G4double GC_PMT_Cone_Min_Inner_Radius=1.25*cm;
	G4double GC_PMT_Cone_Length=3.0185*cm;
	G4double GC_PMT_Cone_Thickness=0.001*cm;


	G4VSolid* GC_PMT_Cone;

	if ( GC_PMT_Cone_Shape=="G4Polyhedra" ) {
		G4int numZPlanes=2;
		G4double* zPlane=new G4double[numZPlanes];
		G4double* rInner=new G4double[numZPlanes];
		G4double* rOuter=new G4double[numZPlanes];
		zPlane[0]=GC_PMT_Cone_Length*0.5;
		zPlane[1]=-GC_PMT_Cone_Length*0.5;
		if ( fabs(GC_PMT_Cone_Length)<1e-5 ) {
			zPlane[0]=1.0; zPlane[1]=-1.0;
		}
		rInner[0]=GC_PMT_Cone_Min_Inner_Radius;
		rInner[1]=GC_PMT_Cone_Max_Inner_Radius;
		rOuter[0]=GC_PMT_Cone_Min_Inner_Radius+GC_PMT_Cone_Thickness;
		rOuter[1]=GC_PMT_Cone_Max_Inner_Radius+GC_PMT_Cone_Thickness;

		GC_PMT_Cone=new G4Polyhedra(GC_PMT_Cone_Name.data(),0,360*deg,6,numZPlanes,zPlane,rInner,rOuter);
	}
	else if ( GC_PMT_Cone_Shape=="G4Cons" ) {
		//GC_PMT_Cone=Construct_Square_Opening_Cone(GC_PMT_Cone_Name,
		//		GC_PMT_Cone_Max_Inner_Radius,GC_PMT_Cone_Max_Inner_Radius+GC_PMT_Cone_Thickness,GC_PMT_Cone_Min_Inner_Radius,GC_PMT_Cone_Min_Inner_Radius+GC_PMT_Cone_Thickness,
		//		GC_PMT_Cone_Length*0.5,0,360*deg,
		//		usermanager->GetGC_PMT_Ver_Space(),usermanager->GetGC_PMT_Hor_Space());
		GC_PMT_Cone=Construct_Square_Opening_Cone(GC_PMT_Cone_Name,
				0,GC_PMT_Cone_Max_Inner_Radius,0,GC_PMT_Cone_Min_Inner_Radius,
				GC_PMT_Cone_Length*0.5,0,360*deg,
				3.1*cm,3.1*cm);
	}

	G4LogicalVolume* GC_PMT_Cone_log = new G4LogicalVolume(GC_PMT_Cone, GetMaterial(G4String("Glass")), GC_PMT_Cone_Name+"_log", 0, 0, 0);
	G4VisAttributes* GC_PMT_Cone_log_VisAtt = new G4VisAttributes();
	GC_PMT_Cone_log_VisAtt->SetColor(GetColor(G4String("Yellow")));
	GC_PMT_Cone_log_VisAtt->SetVisibility(true);
	GC_PMT_Cone_log->SetVisAttributes(GC_PMT_Cone_log_VisAtt);

	// in principle, frame area is (2R*sin(60)*(Num_of_rows-1)+2R)*9*2R and covered area is pi*R^2*((Nint(Num_of_rows*0.5)+Num_of_rows%2)*9+Nint(Num_of_rows*0.5)*8)
	// so active area factor is pi*((Nint(Num_of_rows*0.5)+Num_of_rows%2)*9+Nint(Num_of_rows*0.5)*8)/(36*(sin(60)*(Num_of_rows-1)+1))
	// for Num_of_rows=65, factor=85.53%

	G4AssemblyVolume* assemblyPMT = new G4AssemblyVolume();

	G4int Num_of_rows = 60;
	G4int Num_Of_PMTs_In_Odd_Row= 9;
	G4int Num_Of_PMTs_In_Even_Row= 8;
	G4bool PMT_Is_Sandwiched=true;
	G4int Max_Number_Of_PMTs_In_A_Row=9;
	G4double Space_Between_Rows = 2*GC_PMT_Cover_Radius*sin(60.*deg);
	G4double VSpace=3.1*cm;
	if ( VSpace>Space_Between_Rows ) {
		Space_Between_Rows=VSpace;
	}
	G4double Space_Between_Cols=2*GC_PMT_Cover_Radius;
	G4double HSpace=3.1*cm;
	if ( HSpace>Space_Between_Cols ) {
		Space_Between_Cols=HSpace;
	}

	G4int col;
	G4double xpos,ypos,zpos;
	G4double V_Len=(Num_of_rows-1)*Space_Between_Rows;//vertical length of the center of most top and bottom PMTs
	G4double H_Len=(Max_Number_Of_PMTs_In_A_Row-1)*Space_Between_Cols;//horizontal length of the center of most left and right PMTs

	//in assemblyPMT, the center is the center of PMT
	//xpos=-(GC_PMT_Box_FullSize.x()-GC_PMT_Length)*0.5;//total length of pmt_box is (x gap or gc pmt cone length)+GC_PMT_Length
	//xpos=0;//total length of pmt_box is (x gap or gc pmt cone length)+GC_PMT_Length
	xpos=0;//total length of pmt_box is (x gap or gc pmt cone length)+GC_PMT_Length
	//Translation_PMT(row=0,col=0)=(0,-91.509,-14.3275)
	//Row number is 1-N from bottom
	//Col number is 1-N from left
	for ( G4int row=0; row<Num_of_rows; row++ ) {
		ypos = row*Space_Between_Rows-V_Len*0.5;
		if( row%2 == 0 ) {//odd
			for ( col=0; col<Num_Of_PMTs_In_Odd_Row; col++ ) {
				zpos = (col+0.5*((Num_Of_PMTs_In_Odd_Row-Num_Of_PMTs_In_Odd_Row%2)%2))*Space_Between_Cols-H_Len*0.5;
				rm=rm.IDENTITY;
				rm.rotateY(-90*deg);
				Translation.set(xpos-GC_PMT_Glass_Thickness*0.5,ypos,zpos);
				//printf("row=%d,col=%d,y=%g cm,z=%g cm\n",row+1,col+1,ypos/cm,zpos/cm);
				assemblyPMT->AddPlacedVolume(GC_PMT_log, Translation,&rm);
				Translation.set(xpos,ypos,zpos);
				assemblyPMT->AddPlacedVolume(GC_PMT_Cover_log, Translation,&rm);
				if ( fabs(GC_PMT_Glass_Thickness)>1e-5 ) {
					Translation.set(xpos+(GC_PMT_Length-GC_PMT_Glass_Thickness)*0.5,ypos,zpos);
					assemblyPMT->AddPlacedVolume(GC_PMT_Glass_log, Translation,&rm);
				}
				if ( fabs(GC_PMT_Cone_Length)>1e-5 ) {
					if ( GC_PMT_Cone_Shape=="G4Polyhedra" ) {
						Translation.setX(xpos+GC_PMT_Cone_Length*0.5+GC_PMT_Length*0.5);
						rm.rotateX(30*deg);
					}
					else if ( GC_PMT_Cone_Shape=="G4Cons" ) {
						Translation.setX(xpos+GC_PMT_Cone_Length*0.5+GC_PMT_Length*0.5);
					}
					assemblyPMT->AddPlacedVolume(GC_PMT_Cone_log, Translation,&rm);
				}
			}
		}
		else {//even
			for ( col=0; col<Num_Of_PMTs_In_Even_Row; col++ ) {
				if ( PMT_Is_Sandwiched ) {
					zpos = (col+0.5*((Num_Of_PMTs_In_Even_Row-Num_Of_PMTs_In_Even_Row%2-1)%2))*Space_Between_Cols-H_Len*0.5;
				}
				else {
					zpos = (col+0.5*((Num_Of_PMTs_In_Even_Row-Num_Of_PMTs_In_Even_Row%2)%2))*Space_Between_Cols-H_Len*0.5;
				}
				rm=rm.IDENTITY;
				rm.rotateY(-90*deg);
				Translation.set(xpos-GC_PMT_Glass_Thickness*0.5,ypos,zpos);
				//printf("row=%d,col=%d,y=%g cm,z=%g cm\n",row+1,col+1,ypos/cm,zpos/cm);
				assemblyPMT->AddPlacedVolume(GC_PMT_log, Translation,&rm);
				Translation.set(xpos,ypos,zpos);
				assemblyPMT->AddPlacedVolume(GC_PMT_Cover_log, Translation,&rm);
				if ( fabs(GC_PMT_Glass_Thickness)>1e-5 ) {
					Translation.set(xpos+(GC_PMT_Length-GC_PMT_Glass_Thickness)*0.5,ypos,zpos);
					assemblyPMT->AddPlacedVolume(GC_PMT_Glass_log, Translation,&rm);
				}
				if ( fabs(GC_PMT_Cone_Length)>1e-5 ) {
					if ( GC_PMT_Cone_Shape=="G4Polyhedra" ) {
						Translation.setX(xpos+GC_PMT_Cone_Length*0.5+GC_PMT_Length*0.5);
						rm.rotateX(30*deg);
					}
					else if ( GC_PMT_Cone_Shape=="G4Cons" ) {
						Translation.setX(xpos+GC_PMT_Cone_Length*0.5+GC_PMT_Length*0.5);
					}
					assemblyPMT->AddPlacedVolume(GC_PMT_Cone_log, Translation,&rm);
				}
			}
		}
	}
	//Jail 1-N from bottom to top
	for ( i = 0; i < Num_Of_GC_PMT_Box_Jail_Bars; ++i ) {
		rm=rm.IDENTITY;
		ypos = (-1.0*(Num_Of_GC_PMT_Box_Jail_Bars-1)/2.+i)*Space_Between_Rows;
		Translation.set((GC_PMT_Length+GC_PMT_Box_Jail_FullSize.x())*0.5+GC_PMT_Box_Jail_Space,ypos,0);
		assemblyPMT->AddPlacedVolume(GC_PMT_Box_Jail_log,Translation,&rm);
	}
	rm=rm.IDENTITY;
	Translation.set((GC_PMT_Box_FullSize.x()-GC_PMT_Length)*0.5-GC_PMT_Box_Thickness*0.5,0,0);
	assemblyPMT->AddPlacedVolume(GC_PMT_Box_log, Translation, &rm);

	//rm = rm.IDENTITY;
	//rm.rotateY(usermanager->GetGC_PMT_Box_RotY_Angle());
	//mirror_index=usermanager->GetGC_Num_of_Mirrors()*0.5+1;
	//Translation.set(usermanager->GetGC_Mirror_New_Offset(mirror_index).x(),0,0);
	//newAxis.set(-1,0,0);//-X-axis
	//newAxis.rotateY(usermanager->GetGC_PMT_Box_RotY_Angle());
	//Translation=Translation+newAxis*(usermanager->GetGC_PMT_Box_Distance()+GC_PMT_Length*0.5);


        G4RotationMatrix PACS_rm;
        G4ThreeVector myTranslation,mynewAxis;
        PACS_rm      = PACS_rm.IDENTITY;
        PACS_rm.rotateY(55*deg);
        myTranslation.set(0.0,0,0);
        mynewAxis.set(-1,0,0);//-X-axis
        mynewAxis.rotateY(55*deg);
        myTranslation=myTranslation+mynewAxis*(65*cm+13.75*cm*0.5);//Original center is PMT center, so move PMT_Length*0.5 to surface

        G4ThreeVector PACS_offset=myTranslation;


	

	assemblyPMT->MakeImprint(GC_Tank_log,PACS_offset,&PACS_rm);//Do not use G4VPhysicalVolume otherwise no reflection for pmt_cone

	/* --------------------------------------------------------------------------*/
	//Tank_phys
	Tank_phys = new G4PVPlacement(NULL, G4ThreeVector(), GC_Tank_log, GC_Tank_Name+"_phys", bblog, false, 0);



	G4ThreeVector V1,V2;


        const G4int nEntries=25;
        G4double nm_lambda;
        G4double PhotonEnergy[nEntries];
        G4double PhotonWaveLength[nEntries]={
            650*nm, 600*nm, 550*nm, 500*nm, 450*nm, 400*nm, 350*nm,
            350*nm, 330*nm, 310*nm, 290*nm, 270*nm, 250*nm, 230*nm, 210*nm,
            210*nm, 208*nm, 206*nm, 204*nm, 202*nm, 200*nm,
            199*nm, 195*nm, 190*nm, 185*nm
        };


	/* --------------------------------------------------------------------------*/
	//GC_Mirror surface
	G4double Mirror_Reflectivity[nEntries];
	G4double Mirror_Efficiency[nEntries];
	G4MaterialPropertiesTable** Mirror_SPT=new G4MaterialPropertiesTable*[num_of_mirrors];
	G4OpticalSurface** OpMirrorSurface=new G4OpticalSurface*[num_of_mirrors];
	G4LogicalSkinSurface** GC_Mirror_Surface=new G4LogicalSkinSurface*[num_of_mirrors];

	//hc/(1eV)=1240nm wavelength of 1 eV particle
	//****************************************
	//for 185nm~200nm,E=6.71*eV~6.20*eV
	//p0                        =      57.1302   +/-   4.35232e+07
	//p1                        =     -1.12253   +/-   452877
	//p2                        =   0.00585938   +/-   1176.3
	//
	//****************************************
	//for 200nm~210nm,E=6.199*eV~5.904*eV
	//p0                        =       -13903   +/-   2058.18
	//p1                        =       135.85   +/-   20.0863
	//p2                        =        -0.33   +/-   0.0489898
	//
	//****************************************
	//for 210nm~350nm,E=5.90*eV~3.543*eV
	//p0                        =     -504.252   +/-   43.524
	//p1                        =       5.9471   +/-   0.486514
	//p2                        =   -0.0198695   +/-   0.00178926
	//p3                        =  2.21164e-05   +/-   2.16139e-06
	//
	//****************************************
	//for 350nm~650nm,E=3.54*eV~1.91*eV
	//p0                        =      80.3607   +/-   2.66686
	//p1                        =    0.0553214   +/-   0.0114894
	//p2                        = -6.78571e-05   +/-   1.20585e-05
	for ( i = 0; i < num_of_mirrors; ++i ) {
		mirror_index=i+1;
		for ( j = 0; j < nEntries; ++j ) {
                    nm_lambda=PhotonWaveLength[j]/nm;
                    if ( j<7 ) {
                        Mirror_Reflectivity[j]=80.3607+0.0553214*nm_lambda-6.78571e-05*nm_lambda*nm_lambda;
                    }
                    else if ( j<15 ) {
                        Mirror_Reflectivity[j]=-504.252+5.9471*nm_lambda-0.0198695*nm_lambda*nm_lambda+2.21164e-05*nm_lambda*nm_lambda*nm_lambda;
                    }
                    else if ( j<21 ) {
                        Mirror_Reflectivity[j]=-13903.+135.85*nm_lambda-0.33*nm_lambda*nm_lambda;
                    }
                    else {
                        Mirror_Reflectivity[j]=57.1302-1.12253*nm_lambda+0.00585938*nm_lambda*nm_lambda;
                    }
                    Mirror_Reflectivity[j]/=100.;
                    Mirror_Efficiency[j]=0.;
                    //if ( i==0 ) {
			//	printf("Line %d:lamdba=%g nm, Mirror_Reflectivity[%d]=%g\n",__LINE__,nm_lambda,j,Mirror_Reflectivity[j]);
			//}
		}
		Mirror_SPT[i]=new G4MaterialPropertiesTable();
		Mirror_SPT[i]->AddProperty("REFLECTIVITY", PhotonEnergy, Mirror_Reflectivity, nEntries);
		Mirror_SPT[i]->AddProperty("EFFICIENCY",   PhotonEnergy, Mirror_Efficiency, nEntries);

		sprintf(tmpname,"OpMirrorSurface_%d",mirror_index);
		OpMirrorSurface[i] = new G4OpticalSurface(tmpname);
		OpMirrorSurface[i]->SetType(dielectric_metal);
		OpMirrorSurface[i]->SetFinish(polished);
		OpMirrorSurface[i]->SetModel(unified);
		OpMirrorSurface[i]->SetMaterialPropertiesTable(Mirror_SPT[i]);
		sprintf(tmpname,"GC_Mirror_Surface_%d",mirror_index);
		GC_Mirror_Surface[i]=new G4LogicalSkinSurface(tmpname,GC_Mirror_log[i],OpMirrorSurface[i]);
	}
	/* --------------------------------------------------------------------------*/


	/* --------------------------------------------------------------------------*/
	//GC_PMT surface
	G4double PMT_Reflectivity[nEntries];
	G4double PMT_Efficiency[nEntries];
	for ( i = 0; i < nEntries; ++i ) {
		PMT_Reflectivity[i]=0.;
		PMT_Efficiency[i]=1.0;
	}
	G4MaterialPropertiesTable* PMT_SPT = new G4MaterialPropertiesTable();
	PMT_SPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMT_Reflectivity, nEntries);
	PMT_SPT->AddProperty("EFFICIENCY",   PhotonEnergy, PMT_Efficiency, nEntries);

	G4OpticalSurface* OpPMTSurface = new G4OpticalSurface("OpPMTSurface");
	OpPMTSurface->SetType(dielectric_metal);
	OpPMTSurface->SetFinish(polished);
	OpPMTSurface->SetModel(unified);
	OpPMTSurface->SetMaterialPropertiesTable(PMT_SPT);

	G4LogicalSkinSurface* GC_PMT_Surface;

	GC_PMT_Surface = new G4LogicalSkinSurface("GC_PMT_Surface", GC_PMT_log, OpPMTSurface);

	////GC_PMT_Glass Surface
	//G4OpticalSurface* OpPMTGlassSurface = new G4OpticalSurface("OpPMTGlassSurface");
	//OpPMTGlassSurface->SetType(dielectric_dielectric);
	//OpPMTGlassSurface->SetFinish(polished);
	//OpPMTGlassSurface->SetModel(unified);
	//G4LogicalSkinSurface* GC_PMT_Glass_Surface;

	//GC_PMT_Glass_Surface = new G4LogicalSkinSurface("GC_PMT_Glass_Surface", GC_PMT_Glass_log, OpPMTGlassSurface);

	//GC_PMT_Cover surface
	G4double PMT_Cover_Reflectivity[nEntries];
	G4double PMT_Cover_Efficiency[nEntries];
	for ( i = 0; i < nEntries; ++i ) {
		PMT_Cover_Reflectivity[i]=0.;
		PMT_Cover_Efficiency[i]=0.;
	}
	G4MaterialPropertiesTable* PMT_Cover_SPT = new G4MaterialPropertiesTable();
	PMT_Cover_SPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMT_Cover_Reflectivity, nEntries);
	PMT_Cover_SPT->AddProperty("EFFICIENCY",   PhotonEnergy, PMT_Cover_Efficiency, nEntries);

	G4OpticalSurface* OpPMTCoverSurface = new G4OpticalSurface("OpPMTCoverSurface");
	OpPMTCoverSurface->SetType(dielectric_metal);
	OpPMTCoverSurface->SetFinish(polished);
	OpPMTCoverSurface->SetModel(unified);
	OpPMTCoverSurface->SetMaterialPropertiesTable(PMT_Cover_SPT);

	G4LogicalSkinSurface* GC_PMT_Box_Surface;
	G4LogicalSkinSurface* GC_PMT_Cover_Surface;

	GC_PMT_Box_Surface  = new G4LogicalSkinSurface("GC_PMT_Box_Surface", GC_PMT_Box_log, OpPMTCoverSurface);
	GC_PMT_Cover_Surface = new G4LogicalSkinSurface("GC_PMT_Cover_Surface", GC_PMT_Cover_log, OpPMTCoverSurface);

	//GC_PMT_Cone
	G4double gc_pmt_cone_ref= 0.95;
	G4double PMT_Cone_Reflectivity[nEntries];
	G4double PMT_Cone_Efficiency[nEntries];
	for ( i = 0; i < nEntries; ++i ) {
		PMT_Cone_Reflectivity[i]=gc_pmt_cone_ref;
		PMT_Cone_Efficiency[i]=0.;//if ref<1, then if eff=0, cone will absorb photon, eff=1, cone will det photon
	}
	G4MaterialPropertiesTable* PMT_Cone_SPT = new G4MaterialPropertiesTable();
	PMT_Cone_SPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMT_Cone_Reflectivity, nEntries);
	PMT_Cone_SPT->AddProperty("EFFICIENCY",   PhotonEnergy, PMT_Cone_Efficiency, nEntries);

	G4OpticalSurface* OpPMTConeSurface = new G4OpticalSurface("OpPMTConeSurface");
	OpPMTConeSurface->SetType(dielectric_metal);
	OpPMTConeSurface->SetFinish(polished);
	OpPMTConeSurface->SetModel(unified);
	OpPMTConeSurface->SetMaterialPropertiesTable(PMT_Cone_SPT);

	G4LogicalSkinSurface* GC_PMT_Cone_Surface;

	GC_PMT_Cone_Surface = new G4LogicalSkinSurface("GC_PMT_Cone_Surface", GC_PMT_Cone_log, OpPMTConeSurface);

	//GC_PMT_Box_Jail
	//For dielectric_dielectric
	//Polished: fresnel reflection, total internal reflection and fresnel refraction
	//PolishedFrontPainted: spike reflection and absorption
	//PolishedBackPainted:  spike reflection, lobe reflection, backscatter, lambertian reflection, fresnel refraction and absorption
	//Ground:  spike reflection, lobe reflection, backscatter, lambertian reflection and fresnel refraction
	//GroundFrontPainted: lambertian reflection and absorption
	//GroundBackPainted:  spike reflection, lobe reflection, backscatter, lambertian reflection, fresnel refraction and absorption
	G4double gc_pmt_box_jail_ref=0.8;
	G4double PMT_Box_Jail_Reflectivity[nEntries];
	G4double PMT_Box_Jail_Efficiency[nEntries];
	for ( i = 0; i < nEntries; ++i ) {
		PMT_Box_Jail_Reflectivity[i]=gc_pmt_box_jail_ref;
		PMT_Box_Jail_Efficiency[i]=0;//if ref<1, then if eff=0, cone will absorb photon, eff=1, cone will det photon
	}
	G4MaterialPropertiesTable* PMT_Box_Jail_SPT = new G4MaterialPropertiesTable();
	PMT_Box_Jail_SPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMT_Box_Jail_Reflectivity, nEntries);
	PMT_Box_Jail_SPT->AddProperty("EFFICIENCY",   PhotonEnergy, PMT_Box_Jail_Efficiency, nEntries);

	G4OpticalSurface* OpPMTBoxJailSurface = new G4OpticalSurface("OpPMTBoxJailSurface");
	OpPMTBoxJailSurface->SetType(dielectric_metal);
	OpPMTBoxJailSurface->SetFinish(polished);
	OpPMTBoxJailSurface->SetModel(unified);
	//OpPMTBoxJailSurface->SetPolish(0.8);
	OpPMTBoxJailSurface->SetMaterialPropertiesTable(PMT_Box_Jail_SPT);

	G4LogicalSkinSurface* GC_PMT_Box_Jail_Surface;

	GC_PMT_Box_Jail_Surface = new G4LogicalSkinSurface("GC_PMT_Box_Jail_Surface", GC_PMT_Box_Jail_log, OpPMTBoxJailSurface);

	/* --------------------------------------------------------------------------*/


	return;
}

/* --------------------------------------------------------------------------*/
/**
 * @brief Construct simple shape like G4Box...
 *
 * @param aName Shape Name
 * @param aShape Shape Class Name like G4Box
 * @param aFullSize Full size of dimensions
 *
 * @return G4VSolid
 */
G4VSolid* G4SBSGrinch::ConstructSimple(const G4String& aName, const G4String& aShape, const G4ThreeVector& aFullSize) {
	if ( aShape=="G4Box" ) {
		return new G4Box(aName,aFullSize.x()*0.5,aFullSize.y()*0.5,aFullSize.z()*0.5);
	}
}

/* --------------------------------------------------------------------------*/
/**
 * @brief Get G4Colour for G4VisAttributes
 *
 * @param aColor Color name:black, blue, cyan, gray, green, grey, magenta, red, white, yellow
 *
 * @return color (default is black)
 */
G4Colour G4SBSGrinch::GetColor(const G4String& aColor) {
	G4Colour outColor;
	if ( !aColor.compareTo("black",G4String::ignoreCase) ) {
		outColor=G4Colour::Black();
	}
	else if ( !aColor.compareTo("blue",G4String::ignoreCase) ) {
		outColor=G4Colour::Blue();
	}
	else if ( !aColor.compareTo("cyan",G4String::ignoreCase) ) {
		outColor=G4Colour::Cyan();
	}
	else if ( !aColor.compareTo("gray",G4String::ignoreCase) ) {
		outColor=G4Colour::Gray();
	}
	else if ( !aColor.compareTo("green",G4String::ignoreCase) ) {
		outColor=G4Colour::Green();
	}
	else if ( !aColor.compareTo("grey",G4String::ignoreCase) ) {
		outColor=G4Colour::Grey();
	}
	else if ( !aColor.compareTo("magenta",G4String::ignoreCase) ) {
		outColor=G4Colour::Magenta();
	}
	else if ( !aColor.compareTo("red",G4String::ignoreCase) ) {
		outColor=G4Colour::Red();
	}
	else if ( !aColor.compareTo("white",G4String::ignoreCase) ) {
		outColor=G4Colour::White();
	}
	else if ( !aColor.compareTo("yellow",G4String::ignoreCase) ) {
		outColor=G4Colour::Yellow();
	}
	else {
		outColor=G4Colour::Black();
	}
	return outColor;
}

