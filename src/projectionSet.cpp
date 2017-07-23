#include "projectionSet.h"
#include <iostream>
#include <ostream>
#include <vector>

using namespace std;

		ProjectionSet::ProjectionSet() {}

		ProjectionSet::~ProjectionSet()
		{
			if(mIndexList!=NULL)
				delete[] mIndexList;
		}

		void ProjectionSet::init(unsigned int aProjectionCount)
		{
			mSize = aProjectionCount;
			mIndexList = new unsigned int[mSize];

			for(unsigned int i=0; i<mSize; i++)
				mIndexList[i]=i;

			order();
		}

		void ProjectionSet::init(unsigned int aProjectionCount, std::set<unsigned int> aSkipProjectionList)
		{
			vector<unsigned int> projectionVector;
            projectionVector.reserve(aProjectionCount);

            for (unsigned int i = 0; i < aProjectionCount; i++)
			{
                if(aSkipProjectionList.find(i) == aSkipProjectionList.end())
                {
                    projectionVector.push_back(i);
                }
            }

			if(projectionVector.empty())
            {
				cout << "No projections left to process (out of " << aProjectionCount << ") Check your projection skip list" << endl;
			}

			mSize = projectionVector.size();
            mIndexList = new unsigned int[mSize];
			for(unsigned int i = 0; i < mSize; i++)
            {
				mIndexList[i] = projectionVector[i];
            }

			order();
		}


		unsigned int ProjectionSet::iterator::first()
		{
			return mIndex;
		}

		unsigned int& ProjectionSet::iterator::second()
		{
			return *mList;
		}

		size_t ProjectionSet::getSize()
		{
			return mSize;
		}

		ProjectionSet::iterator ProjectionSet::begin() const
		{
			return iterator(*this);
		}

		ProjectionSet::iterator ProjectionSet::end() const
		{
			return iterator(mSize);
		}

		ProjectionSet::iterator ProjectionSet::rbegin() const
		{
			return iterator(*this, mSize-1);
		}

		ProjectionSet::iterator ProjectionSet::rend() const
		{
			return iterator(-1);
		}

		ProjectionSet::iterator::iterator(const ProjectionSet& aList): mList(aList.mIndexList), mIndex(0) {}

		ProjectionSet::iterator::iterator(const iterator& aIt) : mList(aIt.mList), mIndex(0) {}

		ProjectionSet::iterator::iterator(const iterator& aIt, size_t aSize): mList(&aIt.mList[aSize]), mIndex((unsigned int)aSize) {}

		ProjectionSet::iterator::iterator(size_t aSize): mList(0), mIndex((unsigned int)aSize) {}

		unsigned int& ProjectionSet::iterator::operator*() {return *mList;}

		ProjectionSet::iterator& ProjectionSet::iterator::operator++() {++mList; mIndex++; return *this;}
		ProjectionSet::iterator ProjectionSet::iterator::operator++(int) {iterator tmp(*this); operator++(); return tmp;}

		ProjectionSet::iterator& ProjectionSet::iterator::operator--() {--mList; mIndex--; return *this;}
		ProjectionSet::iterator ProjectionSet::iterator::operator--(int) {iterator tmp(*this); operator--(); return tmp;}

		bool ProjectionSet::iterator::operator==(const iterator& aIt)
		{
			if(mIndex==aIt.mIndex)
				return true;
			else
				return false;
		}

		bool ProjectionSet::iterator::operator!=(const iterator& aIt)
		{
			if(mIndex!=aIt.mIndex)
				return true;
			else
				return false;
		}


		ProjectionSetMaxAngle::ProjectionSetMaxAngle(){}

		ProjectionSetMaxAngle::~ProjectionSetMaxAngle(){}

		void ProjectionSetMaxAngle::order()
		{
			unsigned int* temp = new unsigned int[mSize];

			unsigned int j=(int)(mSize/2);
			unsigned int k=0;

			for(unsigned int i = 0; i < mSize; i++)
			{
				if((i%2)==0)
				{
					temp[i]=mIndexList[j];
					j++;
				}
				else
				{
					temp[i]=mIndexList[k];
					k++;
				}
			}


			for(unsigned int i = 0; i < mSize; i++)
			{
				mIndexList[i]=temp[i];
			}

			delete[] temp;
		}

		ProjectionSetIdentity::ProjectionSetIdentity(){}

		ProjectionSetIdentity::~ProjectionSetIdentity(){}

		void ProjectionSetIdentity::order()
		{

		}
