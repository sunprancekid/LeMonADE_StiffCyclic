#ifndef HistogramGeneralStatistik1D_H
#define HistogramGeneralStatistik1D_H

#include <iostream>

#include <stdint.h>
#include <stdexcept>
#include <vector>
#include <cmath>

#include "StatisticMoment.h"

class HistogramGeneralStatistik1D
{
public:

	// initializer for 1d histogram
	HistogramGeneralStatistik1D()
	{
		nBins = 0;
		nValues= 0.0;
		minVal=0.0;
		maxVal = 0.0;
		binsize = 0.0;
	}

	// contructor for 1d histogram
	HistogramGeneralStatistik1D(double min, double max, int nbins)
	{
		
		nBins = nbins;
		nValues= 0.0;

		minVal=min;
		maxVal = max;
		
		if(nbins>0)
			binsize = (max-min)/double(nbins);
		else
			binsize=0.0;
		
		histogram.resize(nbins);
		histogram_count.resize(nbins);
		
	}

	virtual ~HistogramGeneralStatistik1D()
	{
	}
	
	// resets entire histogram with new minimum, maximum and bin sizes
	void reset(double min=0.0,double max=0.0,size_t nbins=0)
	{
		histogram.clear();

		nBins = nbins;
		nValues= 0.0;

		minVal=min;
		maxVal = max;

		if(nbins>0)
			binsize = (max-min)/double(nbins);
		else
			binsize=0.0;
		
		histogram.resize(nbins);
		histogram = {0.};
		histogram_count.resize(nbins);
		histogram_count = {0};
	}
	
	// geth the bin number (integer) associated with a value
	size_t getBinNo(double x) const
	{
		if(x<minVal)
			throw std::runtime_error("HistogramGeneralStatistik1D::getBinNo()...x value too low\n");
		else if(x>maxVal)
			throw std::runtime_error("HistogramGeneralStatistik1D::getBinNo()...x value too high\n");
		
		return uint32_t( std::floor((x-minVal)/binsize));
	}
	
	// returns the center of a histogram bin that is closest to the value the is passed to the method
	double getClosestBinCenter(double x) const
	{
		return (minVal+getBinNo(x)*binsize + binsize/2.0);
	}
	
	// returns the normalized average associated with bin integer
	double getFirstMomentInBin(size_t bin) const
	{
		return histogram[bin] / nValues;
	}
	
	// double getSecondMomentInBin(size_t bin) const
	// {
	// 	return histogram[bin].ReturnM2();
	// }

	// returns the number of counts associated with a bin
	double getNumCountAtKey(double x) const
	{
		return histogram_count[getBinNo(x)];
	}

	// returns the normalized average associated with double
	double getFirstMomentAtKey(double x) const
	{
		return histogram[getBinNo(x)] / nValues;
	}

	// gets the unnormalized count associated with double
	double getCountAtKey(double x) const
	{
		return histogram[getBinNo(x)];
	}
	
	// returns the value associated with the center of a bin
	double getCenterOfBin(int bin) const
	{
		return (minVal+(bin)*binsize + binsize/2.0);
	}
	
	// returns the value associated with the lower boundary of a bin
	double getBinLowerBoundary(int bin) const
	{
		return (minVal+(bin)*binsize);
	}

	// returns the value associated with the upper boundary of a bin
	double getBinUpperBoundary(int bin) const
	{
		return (minVal+(bin+1)*binsize);
	}

	// return the average moment of the histogram
	double getHistAverage () const
	{
		double avg = 0;
		double sum = 0;
		// loop through each bin, add the statistical contribution of the bins value to the sum
		for (int n = 0; n < nBins; n++) {
			avg += histogram[n] * getCenterOfBin(n);
			sum += histogram_count[n];
		}
		return avg / nValues;
	}

	// return variance of histogram
	double getHistVariance() const
	{
		double avg = getHistAverage();
		double sum = 0;
		double var = 0.;
		// loop through, take the square difference between the average and the bin value
		for (int n = 0; n < nBins; n++) {
			var += histogram[n] * std::pow((getCenterOfBin(n) - avg), 2.);
			sum += histogram_count[n];
		}
		return var / (sum - 1);
	}

	// return standard deviation of histogram
	double getHistSTD () const {return std::sqrt(getHistVariance());}

	// returns the number of bins in histogram
	int getNBins() const
	{
		return nBins;
	}

	// returns the total number of times values have been accumulated within the histogram
	double getNCounts() const
	{
		return nValues;
	}
	
	// returns the width of each histogram bin
	double getBinwidth() const
	{
		return binsize;
	}
	
	// returns the maximum value which can be accepted by the histogram (upper boundary of last bin)
	double getMinCoordinate() const{return minVal;}

	// returns the minimum value which can be accepted by the histogram (lower boundary of first bin)
	double getMaxCoordinate() const{return maxVal;}

	//double getValueSum() const{return valueSum;}

	// increments a histogram bin associated with double (x) by a statistical weight (usually one)
	void addValue(double x, double statisticalWeight=1.0) 
	{
		
		if( ( x > maxVal) || (x<minVal))
		{
			std::stringstream errormessage;
			errormessage<<"HistogramGeneralStatistik1D::addValue(). Value "<<x<<" outside boundaries "<<minVal<<" "<<maxVal<<std::endl;
			throw std::runtime_error(errormessage.str());	
		}
		// increment the histogram bins
		histogram[getBinNo(x)] += statisticalWeight;
		histogram_count[getBinNo(x)] += 1;
		// cummulate total number of values
		nValues+=1.0;
	}
	
	// resets a histogram bin associated with a value (x), including the bin count, to a statistical weight
	void resetValue(double x, double statisticalWeight=1.0)
	{

		if( ( x > maxVal) || (x<minVal))
		{
			std::stringstream errormessage;
			errormessage<<"HistogramGeneralStatistik1D::addValue(). Value "<<x<<" outside boundaries "<<minVal<<" "<<maxVal<<std::endl;
			throw std::runtime_error(errormessage.str());
		}
		// reset the histogram bins, total value count
		histogram[getBinNo(x)] = 0.;
		histogram_count[getBinNo(x)] = 0;
		reset_values_count();
		// accumulate histogram bins
		histogram[getBinNo(x)] += statisticalWeight;
		histogram_count[getBinNo(x)] += 1;
		nValues+=1.0;
	}

	// resets the total number of histogram counts to the count associated with each bin
	void reset_values_count ()
	{
		nValues = 0.;
		for (int n = 0; n < nBins; n++) {
			nValues += (int) histogram_count[nBins];
		}
	}

	// returns the center of each bin in the histogram as an array
	std::vector<double> getVectorBins() const
	{
		std::vector<double> positions;
		positions.resize(nBins);

		for(size_t n=0;n<nBins;n++)
		{
			positions[n]=getCenterOfBin(n);
		}

		return positions;
	}
	
	
	
private:


	size_t nBins; 
	double nValues;
	double minVal; // untere Grenze
	double maxVal; // obere Grenze
	double binsize; //Intervalleinteilung

	std::vector<double> histogram;
	std::vector<int> histogram_count;
};

#endif //HistogramGeneralStatistik1D.h
