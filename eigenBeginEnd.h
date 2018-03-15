#ifndef ___EigenIteratorHelper___
#define ___EigenIteratorHelper___

#if __cplusplus >= 201103L
  /**
   *  @brief  Return an iterator pointing to the first element of
   *          the valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline _Tp*
    begin(Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)); }

  /**
   *  @brief  Return an iterator pointing to the first element of
   *          the const valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline const _Tp*
    begin(const Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)); }

  /**
   *  @brief  Return an iterator pointing to one past the last element of
   *          the valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline _Tp*
    end(Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)) + __va.size(); }

  /**
   *  @brief  Return an iterator pointing to one past the last element of
   *          the const valarray.
   *  @param  __va  valarray.
   */
  template<class _Tp>
    inline const _Tp*
    end(const Eigen::Array<_Tp, -1, 1>& __va)
    { return std::addressof(__va(0)) + __va.size(); }

#endif

#endif
