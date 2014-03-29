/*
parameters
radius of MedianImageFilter = 6
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SimpleITK.h"

#include "itkVector.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageKernelOperator.h"

#include "vnl/vnl_math.h"

#include "itkMedianImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"

#include "itkImageDuplicator.h"

const double PI = 3.14159265359;
namespace sitk = itk::simple;
typedef itk::Image< float, 2 > SliceImageType;


int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " USshoulder.exe"<< std::endl;
    return EXIT_FAILURE;
    }

  // declare typedefs here
  	typedef itk::ImageFileReader<SliceImageType> SliceReader;
	typedef itk::ExtractImageFilter< SliceImageType, SliceImageType > ExtractFilterType;
	typedef itk::MedianImageFilter< SliceImageType, SliceImageType >  MedianFilterType;
	typedef itk::ExtractImageFilter< SliceImageType, SliceImageType > ExtractFilterType;
	typedef itk::AddImageFilter < SliceImageType, SliceImageType > AddFilterType;
	typedef itk::ResampleImageFilter<SliceImageType, SliceImageType> ResampleFilterType;
	typedef itk::TranslationTransform< double, 2 >  TransformType;


	//-----------------------
	// Read and process image
	//-----------------------

	std::string filename = "c:\\miia\\yingyinw\\USshoulder\\shoulderSource\\yy_3_26_2014_36.png";
	std::string filterType = argv[1];
  
	SliceReader::Pointer reader = SliceReader::New();
	reader->SetFileName(filename);
	reader->Update();

	SliceImageType::Pointer inputImage = SliceImageType::New();
	inputImage = reader->GetOutput();
	SliceImageType::SizeType inputImageSize;
	inputImageSize[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
	inputImageSize[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
	
	printf("Original image size %d %d\n", inputImageSize[0], inputImageSize[1]);
	printf("Image spacing %f %f\n", reader->GetOutput()->GetSpacing()[0], reader->GetOutput()->GetSpacing()[1] );
	//sitk::Show( sitk::Image( reader->GetOutput() ) );

	
	//
	// Crop Original Screen Image
	//
	
	// hard code screen borders
	int screenX1 = 192;
	int screenX2 = 588;
	int screenY1 = 31;
	int screenY2 = 411;

	SliceImageType::IndexType screenStart;
	screenStart[0] = screenX1;
	screenStart[1] = screenY1; 
	SliceImageType::SizeType screenSize;
	screenSize[0] = screenX2 - screenX1;
	screenSize[1] = screenY2 - screenY1; 
	SliceImageType::RegionType desiredRegionScreen( screenStart, screenSize );

	ExtractFilterType::Pointer extractFilterScreen = ExtractFilterType::New();
	extractFilterScreen->SetExtractionRegion( desiredRegionScreen );
	extractFilterScreen->SetInput( inputImage );
	extractFilterScreen->SetDirectionCollapseToIdentity(); // This is required.

	//sitk::Show( sitk::Image( extractFilterScreen->GetOutput() ) );
	extractFilterScreen->Update();

	SliceImageType::Pointer croppedImage = SliceImageType::New();
	croppedImage = extractFilterScreen->GetOutput();
	printf("\nExtracted image from screenshot\n");
	printf("Extracted image size %d %d\n", croppedImage->GetLargestPossibleRegion().GetSize()[0], croppedImage->GetLargestPossibleRegion().GetSize()[1]);

	SliceImageType::SizeType croppedImageSize;
	croppedImageSize[0] = croppedImage->GetLargestPossibleRegion().GetSize()[0];
	croppedImageSize[1] = croppedImage->GetLargestPossibleRegion().GetSize()[1];
	//printf("Cropped image size %d %d\n", croppedImageSize[0], croppedImageSize[1]);
	

	typedef itk::ImageDuplicator< SliceImageType > DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(croppedImage);
    duplicator->Update();
    SliceImageType::Pointer clonedImage = duplicator->GetOutput();
	printf("\nCloned image\n");
	printf("Cloned image size %d %d\n", clonedImage->GetLargestPossibleRegion().GetSize()[0], clonedImage->GetLargestPossibleRegion().GetSize()[1]);

	//
	// Crop Top Bright Image
	//
	
	// hard code screen borders
	int topX1 = 115 + screenX1;
	int topX2 = 285 + screenX1;
	int topY1 = 1 + screenY1;
	int topY2 = 181 + screenY1;

	SliceImageType::IndexType topStart;
	topStart[0] = topX1;
	topStart[1] = topY1; 
	SliceImageType::SizeType topSize;
	topSize[0] = topX2 - topX1;
	topSize[1] = topY2 - topY1; 
	SliceImageType::RegionType desiredRegionTop( topStart, topSize );

	ExtractFilterType::Pointer extractFilterTop = ExtractFilterType::New();
	extractFilterTop->SetExtractionRegion( desiredRegionTop );
	extractFilterTop->SetInput( croppedImage );
	extractFilterTop->SetDirectionCollapseToIdentity(); // This is required.

	//sitk::Show( sitk::Image( extractFilterTop->GetOutput() ) );
	extractFilterTop->Update();

	SliceImageType::Pointer topImage = SliceImageType::New();
	topImage = extractFilterTop->GetOutput();
	printf("\nTop image\n");
	printf("Top image size %d %d\n", topImage->GetLargestPossibleRegion().GetSize()[0], topImage->GetLargestPossibleRegion().GetSize()[1]);

	// translate	
	TransformType::Pointer cropShift = TransformType::New();
	TransformType::OutputVectorType translation;
	translation[0] = -int( topX1 );
	translation[1] = -int( topX2 );
	cropShift->Translate( translation );
	//printf("\nBottom section translation x:%f y:%f\n", translation2[0], translation2[1]);		
	
	SliceImageType::PointType newOrigin;
	newOrigin[0] = screenX1+topX1;
	newOrigin[1] = screenY1+topX2;

	// interpolator
	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	ResampleFilterType::Pointer resamplerTop = ResampleFilterType::New();
	resamplerTop->SetTransform( cropShift );	
	resamplerTop->SetInterpolator( interpolator );	
	resamplerTop->SetSize( croppedImageSize );
	resamplerTop->SetOutputOrigin( newOrigin );
	resamplerTop->SetDefaultPixelValue( 0 );
	resamplerTop->SetInput( extractFilterTop->GetOutput() );
	
	//sitk::Show( sitk::Image( resampler->GetOutput() ) );
	resamplerTop->Update();	

	SliceImageType::Pointer topImageShifted = resamplerTop->GetOutput();

	//
	// Crop Bottom Bright Image
	//
	
	// hard code screen borders
	int botX1 = 98 + screenX1;
	int botX2 = 303 + screenX1;
	int botY1 = 181 + screenY1;
	int botY2 = 380 + screenY1;

	SliceImageType::IndexType botStart;
	botStart[0] = botX1;
	botStart[1] = botY1; 
	SliceImageType::SizeType botSize;
	botSize[0] = botX2 - botX1;
	botSize[1] = botY2 - botY1; 
	SliceImageType::RegionType desiredRegionBot( botStart, botSize );

	ExtractFilterType::Pointer extractFilterBot = ExtractFilterType::New();
	extractFilterBot->SetExtractionRegion( desiredRegionBot );
	extractFilterBot->SetInput( croppedImage );
	extractFilterBot->SetDirectionCollapseToIdentity(); // This is required.

	//sitk::Show( sitk::Image( extractFilterBot->GetOutput() ) );
	extractFilterBot->Update();

	SliceImageType::Pointer botImage = SliceImageType::New();
	botImage = extractFilterBot->GetOutput();
	printf("\nBottom image\n");
	printf("Bottom image size %d %d\n", botImage->GetLargestPossibleRegion().GetSize()[0], botImage->GetLargestPossibleRegion().GetSize()[1]);

	ResampleFilterType::Pointer resamplerBot = ResampleFilterType::New();
	resamplerBot->SetTransform( cropShift );	
	resamplerBot->SetInterpolator( interpolator );	
	resamplerBot->SetSize( croppedImageSize );
	resamplerBot->SetOutputOrigin( newOrigin );
	resamplerBot->SetDefaultPixelValue( 0 );
	resamplerBot->SetInput( extractFilterBot->GetOutput() );
	
	//sitk::Show( sitk::Image( resampler->GetOutput() ) );
	resamplerBot->Update();	

	SliceImageType::Pointer botImageShifted = resamplerBot->GetOutput();

	//
	// Create Middle Section
	//
	AddFilterType::Pointer addFilter = AddFilterType::New();
	addFilter->SetInput1( topImageShifted );
	addFilter->SetInput2( botImageShifted );
	addFilter->Update();

	SliceImageType::Pointer brightImage = SliceImageType::New();
	brightImage = addFilter->GetOutput();
	//sitk::Show( sitk::Image( brightImage ) );
	printf("\nTop and bottom merged\n");

	//
	// Set crop region to black 0
	//
	
	SliceImageType::IndexType topStart2;
	topStart2[0] = topX1;//+screenX1;
	topStart2[1] = topY1;//+screenY1; 
	SliceImageType::SizeType topSize2;
	topSize2[0] = topX2-topX1;
	topSize2[1] = topY2-topY1;
 	SliceImageType::RegionType topRegion;
	topRegion.SetSize( topSize2 );
	topRegion.SetIndex( topStart2 );
	
	
	itk::ImageRegionIterator<SliceImageType> topImageIterator( clonedImage,topRegion );
 
	while(!topImageIterator.IsAtEnd())
    {
		topImageIterator.Set( 0 ); 
		//SliceImageType::PixelType val = topImageIterator.Get();
		//std::cout << (int)val << std::endl;
		++topImageIterator;
	}

	printf("\nDone iterating\n");
	printf("\nDark image\n");
	printf("Dark image size %d %d\n", clonedImage->GetLargestPossibleRegion().GetSize()[0], clonedImage->GetLargestPossibleRegion().GetSize()[1]);

	sitk::Image sitkImageOut = sitk::Image( clonedImage );
	sitk::Show( sitkImageOut );
	

	/*
	SliceImageType::IndexType botStart;
	botStart[0] = botX1+newOrigin[0];
	botStart[1] = botY1+newOrigin[1];
	SliceImageType::SizeType botSize;
 	botSize[0] = botX2-botX1;
	botSize[1] = botY2-botY1;
 	SliceImageType::RegionType botRegion;
	botRegion.SetSize( botSize );
	botRegion.SetIndex( botStart );
	
	itk::ImageRegionIterator<SliceImageType> botImageIterator( clonedImage,desiredRegionBot ); //botRegion );
 
	while(!botImageIterator.IsAtEnd())
    {
		SliceImageType::PixelType pixel = 0;		
		botImageIterator.Set( pixel ); 
		++botImageIterator;
	}

	sitk::Show( sitk::Image( clonedImage ) );	

	printf("\nSides cropped out.\n");

	*/
	/*
	//
	// MedianImageFilter 
	//
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();  
	medianFilter->SetInput( inputImage ); 
	medianFilter->SetRadius( atoi( argv[2] ) );
	medianFilter->Update();
  
  
	SliceImageType::Pointer medianImage = medianFilter->GetOutput();
	printf("medianImage size %d %d\n", medianImage->GetLargestPossibleRegion().GetSize()[0], medianImage->GetLargestPossibleRegion().GetSize()[1] );
	printf("medianImage spacing %f %f\n", medianImage->GetSpacing()[0], medianImage->GetSpacing()[1] );
	printf("medianImage origin %f %f\n", medianImage->GetOrigin()[0], medianImage->GetOrigin()[1] );
   
	// Show image
	sitk::Show( sitk::Image( medianFilter->GetOutput() ) );
	*/

  return EXIT_SUCCESS;
}






