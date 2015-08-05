/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hdf5SurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"

#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(hdf5SurfaceWriter);
    //addToRunTimeSelectionTable(surfaceWriter, ensightSurfaceWriter, wordDict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

    template<>
    void Foam::hdf5SurfaceWriter::writeData
    (
        const word& surfaceName,
        const word& fieldName,
        const word& time,
        const Field<scalar>& values
    )
    {

        ioScalar* scalarData;
        scalarData = new ioScalar[values.size()];

        // Loop through the field and construct the array
        forAll(values, iter)
        {
            scalarData[iter] = values[iter];
        }

    char hdffileName[80];
    sprintf
    (
        hdffileName,
        "%s.h5",
        surfaceName.c_str()
    );

    hid_t file_id = H5Fopen(hdffileName, H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t group = H5Gcreate2(file_id, time.c_str(), H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    hsize_t dimsf[1];
    dimsf[0] = values.size();
    hid_t dataspace = H5Screate_simple(1, dimsf, NULL);

    hid_t datatype = H5Tcopy(H5T_SCALAR);

    char datasetName[80];
    sprintf
    (
        datasetName,
        "%s/%s",
        time.c_str(),
        fieldName.c_str()
    );

    hid_t  dataset = H5Dcreate2(file_id, datasetName, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset, H5T_SCALAR, H5S_ALL, 
          H5S_ALL, H5P_DEFAULT, scalarData);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Gclose(group);
    H5Fclose(file_id);

    delete [] scalarData;

    }


    template<>
    void Foam::hdf5SurfaceWriter::writeData
    (
        const word& surfaceName,
        const word& fieldName,
        const word& time,
        const Field<vector>& values
    )
    {
        ioScalar* vectorData;
        vectorData = new ioScalar[values.size()*3];

        // Loop through the field and construct the array
        forAll(values, iter)
        {
            vectorData[3*iter+0] = values[iter].x();
            vectorData[3*iter+1] = values[iter].y();
            vectorData[3*iter+2] = values[iter].z();
        }

    char hdffileName[80];
    sprintf
    (
        hdffileName,
        "%s.h5",
        surfaceName.c_str()
    );

    hid_t file_id = H5Fopen(hdffileName, H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t group = H5Gcreate2(file_id, time.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimsf[2];
    dimsf[0] = values.size();
    dimsf[1] = 3;
    hid_t dataspace = H5Screate_simple(2, dimsf, NULL);

    hid_t datatype = H5Tcopy(H5T_SCALAR);

    char datasetName[80];
    sprintf
    (
        datasetName,
        "%s/%s",
        time.c_str(),
        fieldName.c_str()
    );

    hid_t  dataset = H5Dcreate2(file_id, datasetName, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(dataset, H5T_SCALAR, H5S_ALL, 
          H5S_ALL, H5P_DEFAULT, vectorData);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Gclose(group);
    H5Fclose(file_id);

    delete [] vectorData;
    }


    template<>
    void Foam::hdf5SurfaceWriter::writeData
    (
        const word& surfaceName,
        const word& fieldName,
        const word& time,
        const Field<sphericalTensor>& values
    )
    {
        /*os  << "1 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const sphericalTensor& v = values[elemI];
            os  << float(v[0]) << nl;
        }*/
    }


    template<>
    void Foam::hdf5SurfaceWriter::writeData
    (
        const word& surfaceName,
        const word& fieldName,
        const word& time,
        const Field<symmTensor>& values
    )
    {
        /*os  << "6 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const symmTensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << nl;

        }*/
    }


    template<>
    void Foam::hdf5SurfaceWriter::writeData
    (
        const word& surfaceName,
        const word& fieldName,
        const word& time,
        const Field<tensor>& values
    )
    {
        /*os  << "9 " << values.size() << " float" << nl;

        forAll(values, elemI)
        {
            const tensor& v = values[elemI];
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << ' '
                << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
                << nl;
        }*/
    }

}


// Write generic field in vtk format
template<class Type>
void Foam::hdf5SurfaceWriter::writeData
(
    const word& surfaceName,
    const word& fieldName,
    const word& time,
    const Field<Type>& values
)
{
    /*os  << "1 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        os  << float(0) << nl;
    }*/
}

template<class Type>
void Foam::hdf5SurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << surfaceDir << endl;
    }

    wordList dirComp=outputDir.components();
    word time=dirComp[dirComp.size()-1];

    // geometry should already have been written
    // Values to separate directory (e.g. "scalarField/p")

    fileName foamName(pTraits<Type>::typeName);
    fileName valuesDir(surfaceDir  / (foamName + Field<Type>::typeName));

    if (!isDir(valuesDir))
    {
        mkDir(valuesDir);
    }

    
    // Write data
    OFstream os(outputDir/fieldName + '_' + surfaceName + ".vtk");
    writeData(surfaceName,fieldName,time, values);
    //OFstream(valuesDir/fieldName)()  << values;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hdf5SurfaceWriter::hdf5SurfaceWriter()
:
    surfaceWriter()
{
        // Find and create filename

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hdf5SurfaceWriter::~hdf5SurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hdf5SurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }


    if (verbose)
    {
    Info<< "Not writing geometry to " << surfaceDir << endl;
    }

    char hdffileName[80];
    sprintf
    (
        hdffileName,
        "%s.h5",
        surfaceName.c_str()
    );

    if (!isFile(hdffileName))
    {
    Info<< "Creating new H5 file " << hdffileName << endl;
    hid_t file_id = H5Fcreate(hdffileName, H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose(file_id); 
    }

    hid_t file_id = H5Fopen(hdffileName, H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t status = H5Eset_auto2(H5E_DEFAULT, NULL,NULL);
    status = H5Gget_objinfo (file_id, "/mesh", 0, NULL);
    if (status != 0)
    {
        Info<< "writing mesh to H5 " << hdffileName << endl;
        hid_t group = H5Gcreate2(file_id, "/mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ioScalar* vectorData;
        vectorData = new ioScalar[points.size()*3];

        // Loop through the field and construct the array
        forAll(points, iter)
        {
            vectorData[3*iter+0] = points[iter].x();
            vectorData[3*iter+1] = points[iter].y();
            vectorData[3*iter+2] = points[iter].z();
        }

        hsize_t dimsf[2];
        dimsf[0] = points.size();
        dimsf[1] = 3;
        hid_t dataspace = H5Screate_simple(2, dimsf, NULL);

        hid_t datatype = H5Tcopy(H5T_SCALAR);
        hid_t dataset = H5Dcreate2(file_id, "/mesh/points", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset, H5T_SCALAR, H5S_ALL, 
              H5S_ALL, H5P_DEFAULT, vectorData);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);
 

        delete [] vectorData;

        vectorData = new ioScalar[faces.size()*3];

        // Loop through the field and construct the array
        forAll(faces, iter)
        {
            vectorData[3*iter+0] = faces[iter][0];
            vectorData[3*iter+1] = faces[iter][1];
            vectorData[3*iter+2] = faces[iter][2];
        }


        dimsf[0] = faces.size();
        dimsf[1] = 3;
        dataspace = H5Screate_simple(2, dimsf, NULL);

        datatype = H5Tcopy(H5T_SCALAR);
        dataset = H5Dcreate2(file_id, "/mesh/faces", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset, H5T_SCALAR, H5S_ALL, 
              H5S_ALL, H5P_DEFAULT, vectorData);
        H5Dclose(dataset);
        H5Sclose(dataspace);
        H5Tclose(datatype);


        delete [] vectorData;


        H5Gclose(group);
    }

    status = H5Eset_auto2(H5E_DEFAULT, NULL,NULL);
    H5Fclose(file_id);

    /*
    // Points
    OFstream(surfaceDir/"points")() << points;

    // Faces
    OFstream(surfaceDir/"faces")() << faces;

    // Face centers. Not really necessary but very handy when reusing as inputs
    // for e.g. timeVaryingMapped bc.
    pointField faceCentres(faces.size(),point::zero);

    forAll(faces, faceI)
    {
        faceCentres[faceI] = faces[faceI].centre(points);
    }

    OFstream(surfaceDir/"faceCentres")() << faceCentres;
*/

}


// create write methods
defineSurfaceWriterWriteFields(Foam::hdf5SurfaceWriter);


// ************************************************************************* //
