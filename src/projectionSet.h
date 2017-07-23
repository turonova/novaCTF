#pragma once

#include <set>
#include "common.h"


	class ProjectionSet
	{

    protected:

		unsigned int* mIndexList;	//ordered list of projection indices
		size_t mSize;//number of elements in the projection index list

    public:

		ProjectionSet();

		virtual ~ProjectionSet() = 0;

		void init(unsigned int aProjectionCount);

		void init(unsigned int aProjectionCount, std::set<unsigned int> aSkipProjectionList);

	private:

		virtual void order() = 0;

	public:
		class iterator;
		friend class iterator;

		class iterator{
			unsigned int* mList;
			unsigned int mIndex;		//index of currently used projection

		public:

			iterator(const ProjectionSet& aSet);

			iterator(const iterator& aIt);

			iterator(const iterator& aIt, size_t aSize);

			iterator(size_t aSize);

			unsigned int& operator*();

			iterator& operator++();
			iterator operator++(int);

			iterator& operator--();
			iterator operator--(int);

			bool operator==(const iterator& aIt);

			bool operator!=(const iterator& aIt);

			unsigned int first();

			unsigned int& second();

		};

		size_t getSize();

		iterator begin() const;

		iterator end() const;

		iterator rbegin() const;

		iterator rend() const;

	};


	class ProjectionSetMaxAngle: public ProjectionSet
	{

	public:

		ProjectionSetMaxAngle();

		~ProjectionSetMaxAngle();

	private:

		void order();

	};

	class ProjectionSetIdentity: public ProjectionSet
	{

	public:

		ProjectionSetIdentity();

		~ProjectionSetIdentity();

	private:

		void order();

	};
