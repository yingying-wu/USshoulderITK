/*!
 * \file   reconstruct.cpp
 * \author Ved Vyas <ved@cmu.edu>
 * \date   Mar 2014
 *
 * \brief  Reconstruction of 3D volumetric image from 2D slices with poses
 */

// System includes
#include <algorithm>

// External includes
#include "itkEuler3DTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
//#include "itkCropImageFilter.h"
//#include "itkMedianImageFilter.h"

// Project includes
#include "itkTransformForwardingInterpolateImageFunction.h"

// NOTE: still just demo code but demonstrates most of the workflow
int main(const int argc, const char *argv[])
{

  if( argc < 5 )
  {
  std::cerr << "Usage: " << std::endl;
  std::cerr << argv[0] << "reconstruct.exe \
						  inputfilepath \
						  outputfilepath \
						  lowerCropBoundary \
						  upperCropBoundary"<< std::endl;
  return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;
  typedef itk::Image<PixelType, Dimension>  SliceImageType;

  typedef itk::ImageFileReader<> SliceReader;
  SliceReader::Pointer reader = SliceReader::New();

  printf("Going to load: %s\n", argv[1]);

  //if (argc < 2) {return 1;}

  reader->SetFileName(argv[1]);
  reader->Update();

  /*
  SliceImageType::SizeType cropSize;
  cropSize[0] = atoi(argv[3]);
  cropSize[1] = atoi(argv[4]);

  typedef itk::CropImageFilter <SliceImageType, SliceImageType> CropImageFilterType;
  CropImageFilterType::Pointer cropFilter = CropImageFilterType::New();

  cropFilter->SetInput(reader.GetOutput());
  cropFilter->SetBoundaryCropSize(cropSize);

  cropFilter->Update();
  
  typedef itk::MedianImageFilter< SliceImageType, SliceImageType> MedianFilterType;
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();  
  medianFilter->SetInput( cropFilter->GetOutput() ); 
  medianFilter->SetRadius( 12 );
  medianFilter->Update();
  
  GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();

  gaussianFilter->SetSigma( 0.04 ); 
  gaussianFilter->SetDirection( 0 ); // (theta, r)
  gaussianFilter->SetNormalizeAcrossScale( true );

  gaussianFilter->SetInput( polarResampleFilter->GetOutput() );
  gaussianFilter->Update();
  */
        
  typedef itk::Euler3DTransform<double> TransformType;
  TransformType::Pointer xfm = TransformType::New();

  // xfm->SetCenter(...->GetCenter());
  xfm->SetIdentity();
  xfm->SetRotation(-1.57/2., 0., 0.); // 0.02, 0.03);

  typedef itk::ResampleImageFilter<SliceImageType, SliceImageType>
    ResampleFilterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetTransform(xfm);
  resampler->SetInput(reader->GetOutput());
  //resampler->SetInput(cropFilter->GetOutput());

  typedef itk::LinearInterpolateImageFunction<SliceImageType, double>
    InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  // typedef itk::TransformForwardingInterpolateImageFunction<
  //   itk::NearestNeighborInterpolateImageFunction, SliceImageType, double>
  //   XFMInterpolatorType;
  // XFMInterpolatorType::Pointer xfminterp = XFMInterpolatorType::New();

  resampler->SetInterpolator(interpolator);
  // resampler->SetInterpolator(xfminterp);
  resampler->SetDefaultPixelValue(0);

  const SliceImageType::Pointer inputImage = reader->GetOutput();
  //const SliceImageType::Pointer inputImage = cropFilter->GetOutput();

  const SliceImageType::SpacingType &spacing = inputImage->GetSpacing();
  const SliceImageType::PointType &origin = inputImage->GetOrigin();
  SliceImageType::SizeType size =
    inputImage->GetLargestPossibleRegion().GetSize();

  printf("Size: %d %d %d\n", size[0], size[1], size[2]);
  size[0] *= 2;
  size[1] *= 2;
  size[2] = std::max(size[0], size[1]);

  resampler->SetOutputOrigin(origin);
  resampler->SetOutputSpacing(spacing);
  resampler->SetOutputDirection(inputImage->GetDirection());
  resampler->SetSize(size);

  resampler->Update();

  typedef itk::ImageFileWriter<SliceImageType> VolumeWriter;
  VolumeWriter::Pointer writer = VolumeWriter::New();

  printf("Going to write: %s\n", argv[2]);

  if (argc < 3) {return 2;}

  writer->SetFileName(argv[2]);
  writer->SetInput(resampler->GetOutput());

  writer->Update();

  printf("Done\n");

  return 0;
}
