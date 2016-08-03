// C++ part of the adaptor for paraview catalyst for b2.5 simulation.
// Author: Jure Bartol
// Created on: 22.07.2016
// Modified on: 03.08.2016

#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkImageData.h"
#include "vtkCPPythonAdaptorAPI.h"

#include <vtkPoints.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>


extern "C" void creategrid_(double* crx, double* cry, int* ncrx, int* nx, int* ny)
//crx - x coordinate, cry - y coordinate, ncrx - number of points,
//nx - number of cells in x direction, ny - number of cells in y direction
{
  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
    {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
    }


  //create vtk points
  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(*ncrx);
  int numOfCells = (*nx+2)*(*ny+2);
  for (vtkIdType id = 0; id < numOfCells-1 ; ++id){
    points->SetPoint(0+id*4, crx[id], cry[id], 0.0);
    points->SetPoint(1+id*4, crx[id+numOfCells], cry[id+numOfCells], 0.0);
    points->SetPoint(2+id*4, crx[id+2*numOfCells], cry[id+2*numOfCells], 0.0);
    points->SetPoint(3+id*4, crx[id+3*numOfCells], cry[id+3*numOfCells], 0.0);
}

  //create vtk cells
  vtkPolyData* grid = vtkPolyData::New();
  vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input")->SetGrid(grid);
  grid->SetPoints(points);
  points->Delete();
  grid->Allocate(numOfCells);
  vtkIdType ids[4];
  for (int i = 0 ; i < numOfCells ; i++){
    ids[0] = 0+i*4;
    ids[1] = 1+i*4;
    ids[2] = 3+i*4;
    ids[3] = 2+i*4;
    grid->InsertNextCell(9,4,ids);
  }
}

extern "C" void adddata_(double* scalars, char* name, int* numCells)
{
  vtkCPInputDataDescription* idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->GetInputDescriptionByName("input");

  vtkPolyData* grid = vtkPolyData::SafeDownCast(idd->GetGrid());

  if (!grid)
    {
    vtkGenericWarningMacro("No adaptor grid to attach field data to.");
    return;
    }

  vtkDoubleArray* cellDataDoubleArray = vtkDoubleArray::New();
  cellDataDoubleArray->SetName(name);
  cellDataDoubleArray->SetNumberOfTuples(*numCells);
  cellDataDoubleArray->SetNumberOfComponents(1);
  for (vtkIdType i = 0; i < *numCells; ++i){ 
    cellDataDoubleArray->InsertTuple1(i,scalars[i]);
  }
  grid->GetCellData()->AddArray(cellDataDoubleArray);
  cellDataDoubleArray->Delete();
}
