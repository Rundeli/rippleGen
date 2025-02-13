
import sympy as sp

def setProperties(entry,value,map):
    
    if entry in map:
        if entry in ["xlimit","ylimit"]:
            map[entry] = [float(v) for v in value]
        elif entry in ["offset","particle_diameter","symmetry_axis"]:
            map[entry] = float(value[0])
        elif entry in ["is_period_function","is_symmetry","is_piecewise_function"]:
            if value[0].lower() == 'true':
                map[entry] = True
            elif value[0].lower() == 'false':
                map[entry] = False
        elif entry in ["surface_function_z"]:
            map[entry] = value

def parseSurfaceFunction(functionInString):
    x = sp.Symbol('x')
    parseToSmpy=sp.sympify(functionInString)
    f = sp.lambdify(x,parseToSmpy,'numpy')
    return f

def readInputFile(fileName,config):
    with open(fileName,"r",encoding="utf-8") as f_input:
        for line in f_input:
            if line.startswith("#"):continue
            if not line:continue

            parts = line.split()
            if(len(parts)>=2):
                key = parts[0]
                value = parts[1:]
                setProperties(key,value,config)

def calculateParticlePosition(nParticle_x,nParticle_y,offset,surfaceFunction,particle_diameter):

    particle_Cords = []

    for i in range(nParticle_x):
            
        xi=offset + particle_diameter/2+i*(offset + particle_diameter)

        for m in range(nParticle_y):

            ym=offset + particle_diameter/2+m*(offset + particle_diameter)
            nParticle_z=int(round(surfaceFunction(xi)[0]/particle_diameter))+1

            for n in range(nParticle_z):

                zn=0.5*particle_diameter+n*(offset+particle_diameter)
                particle_Cords.append(f"({xi:.6f} {ym:.6f} {zn:.6f})")

    return particle_Cords

def calculateSymmParticlePosition(nParticle_x,nParticle_y,offset,surfaceFunction,particle_diameter,symmetryAxis):

    particle_Cords=[]

    xi0=symmetryAxis-particle_diameter-offset

    for i in range(nParticle_x):

        xi=xi0-i*(particle_diameter+offset)

        symmetry_x = 2 * symmetryAxis - xi

        for m in range(nParticle_y):

            ym=offset + particle_diameter/2+m*(offset + particle_diameter)
            nParticle_z=int(round(surfaceFunction(xi)[0]/particle_diameter))+1

            for n in range(nParticle_z):
                    
                zn=0.5*particle_diameter+n*(offset+particle_diameter)
                particle_Cords.append(f"({xi:.6f} {ym:.6f} {zn:.6f})")
                particle_Cords.append(f"({symmetry_x:.6f} {ym:.6f} {zn:.6f})")

    #for coord in particle_Cords:
    
        #coords = coord.strip("()").split()
        #x, y, z = float(coords[0]), float(coords[1]), float(coords[2])
        
        #symmetry_x = 2 * symmetryAxis - x
        #particle_Cords.append(f"({symmetry_x:.6f} {y:.6f} {z:.6f})")

    nParticle_zCenter=int(round(surfaceFunction(symmetryAxis)[0]/particle_diameter))+1

    for m in range(nParticle_y):

        ym=offset + particle_diameter/2+m*(offset + particle_diameter)

        for n in range(nParticle_zCenter):
                    
            zn=0.5*particle_diameter+n*(offset+particle_diameter)
            particle_Cords.append(f"({symmetryAxis:.6f} {ym:.6f} {zn:.6f})")
    
    return particle_Cords

def periodShift(particle_Cords,period):

    #not available yet

    return particle_Cords

    

def generateParticles(config):

    xLimit=config["xlimit"]
    yLimit=config["ylimit"]
    offset=config["offset"]
    particle_diameter=config["particle_diameter"]
    isSym=config["is_symmetry"]
    symmetry_axis=config["symmetry_axis"]
    isPeriod=config["is_period_function"]
    period=config["period"]

    surfaceFunction = parseSurfaceFunction(config["surface_function_z"])

    if isSym:

        nParticle_x = int((symmetry_axis-xLimit[0]-0.5*particle_diameter-offset)/(particle_diameter+offset))
        nParticle_y = int((yLimit[1]-yLimit[0])/(particle_diameter+offset))

        particle_Cords=calculateSymmParticlePosition(nParticle_x,nParticle_y,offset,surfaceFunction,particle_diameter,symmetry_axis)

    else:
        nParticle_x = int((xLimit[1]-xLimit[0]-offset)/(particle_diameter+offset))
        nParticle_y = int((yLimit[1]-yLimit[0])/(particle_diameter+offset))

        particle_Cords=calculateParticlePosition(nParticle_x,nParticle_y,offset,surfaceFunction,particle_diameter)

    if isPeriod:

        particle_Cords=periodShift(particle_Cords,period)

    return particle_Cords


def writeOutputFile(particle_Cords):

    with open("kinematicCloudPositions","w") as f_output:
        header = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  | Generated with rippleGen by Rundeli 2025        |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    object      kinematicCloudPositions;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

(
"""
        footer = """\
)

// ************************************************************************* //
"""
        f_output.write(header)
        f_output.write("\n".join(particle_Cords))
        f_output.write("\n".join(footer))

    print("粒子位置文件已成功生成!")
    print("共有"+ f"{len(particle_Cords)}" +"个粒子.")
    return

def main():

    xLimit_=[]
    yLimit_=[]
    offset_=0.0
    is_period_function_=False
    period_=0.0
    is_symmetry_=False
    symmetry_axis_=0.0
    is_piecewise_function_=False
    surface_function_z_=None
    particleDiameter_=0.0

    config={
        "xlimit":xLimit_,
        "ylimit":yLimit_,
        "offset":offset_,
        "is_period_function":is_period_function_,
        "period":period_,
        "is_symmetry":is_symmetry_,
        "symmetry_axis":symmetry_axis_,
        "is_piecewise_function":is_piecewise_function_,
        "surface_function_z":surface_function_z_,
        "particle_diameter":particleDiameter_
    }

    readInputFile("input",config)

    particleCords_= generateParticles(config)

    writeOutputFile(particleCords_)

if __name__ == "__main__":
    main()




    

    

    




   











