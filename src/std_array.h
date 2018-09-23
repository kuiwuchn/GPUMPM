#ifndef __std_array_h__
#define __std_array_h__

#include <cassert>
#include <SPGrid_Utilities.h>

namespace SPGrid{

template<class T,int d>
class std_array{

    template<class T1> struct Make_Pointer_To_Const;
    template<class T1> struct Make_Pointer_To_Const<const T1*> { typedef const T1* type; };
    template<class T1> struct Make_Pointer_To_Const<T1*> { typedef const T1* type; };

public:

    T data[d];
    
    std_array()
    {for(int i=0;i<d;i++) data[i]=T();}
    
    std_array(T a)
    {for(int i=0;i<d;i++) data[i]=a;}

    std_array(T a, T b)
    {
        Static_Assert(d==2);
        data[0]=a;
        data[1]=b;
    }
    std_array(T a, T b, T c)
    {
        Static_Assert(d==3);
        data[0]=a;
        data[1]=b;
        data[2]=c;
    }
    
    std_array(const std_array<int,d>& arr)
    {
        for(int i=0;i<d;i++)
            data[i]=T(arr.data[i]);
    }
        

    template<class T_ARRAY>
    explicit std_array(const T_ARRAY& array)
    {typedef typename Make_Pointer_To_Const<typename T_ARRAY::iterator>::type const_iterator;
    int i=0;for(const_iterator it=array.begin();it!=array.end();it++,i++) data[i]=*it;}

    template<class T_ARRAY>
    T_ARRAY Cast()
    {typedef typename T_ARRAY::iterator iterator;
    T_ARRAY array;int i=0;for(iterator it=array.begin();it!=array.end();it++,i++) *it=data[i];
    return array;}

    std_array operator+(const T& value) const 
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]+value;return result;}

    std_array operator-(const T& value) const 
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]-value;return result;}
    
    std_array operator%(const T& value) const 
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]%value;return result;}

    std_array operator>>(const int value) const 
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]>>value;return result;}

    std_array operator<<(const int value) const 
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]<<value;return result;}

    bool operator==(const std_array& array) const
    {for(int i=0;i<d;i++) if(array.data[i]!=data[i]) return false; return true;}

    bool operator!=(const std_array& array) const
    {for(int i=0;i<d;i++) if(array.data[i]!=data[i]) return true; return false;}

    std_array operator%(const std_array& array) const
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]%array.data[i];return result;}
    
    std_array operator+(const std_array& array) const
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]+array.data[i];return result;}
    
    std_array operator-(const std_array& array) const
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]-array.data[i];return result;}
    
    std_array operator*(const std_array& array) const
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]*array.data[i];return result;}
    
    std_array operator/(const std_array& array) const
    {std_array result;for(int i=0;i<d;i++) result.data[i]=data[i]/array.data[i];return result;}

    T& operator()(const int i)
    {assert(i<d);return data[i];}

    T operator()(const int i) const
    {assert(i<d);return data[i];}

};

template<class T,int d>
std::ostream& operator<<(std::ostream& out,const std_array<T,d>& arr)
{
    out << "[ ";
    for (int i=0;i<d;i++)
        out<<arr.data[i]<<" ";
    out << "]";
    return out;
}

}
#endif
