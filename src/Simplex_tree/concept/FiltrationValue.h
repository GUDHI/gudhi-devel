/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/** \brief Value type for a filtration function on a cell complex.
  *
  * A <EM>filtration</EM> of a cell complex (see FilteredComplex) is
  * a function \f$f:\mathbf{K} \rightarrow \mathbb{R}\f$ satisfying \f$f(\tau)\leq 
  * f(\sigma)\f$ whenever \f$\tau \subseteq \sigma\f$. Ordering the simplices 
  * by increasing filtration values (breaking ties so as a simplex appears after 
  * its subsimplices of same filtration value) provides an indexing scheme 
  * (see IndexingTag).
  */
  struct FiltrationValue {
    /**
     * @brief Has to construct the default value of FiltrationValue. E.g., 0 for a numerical value or {} for a vector.
     */
    FiltrationValue(0);

    /** \brief Operator < is a StrictWeakOrdering. */
    bool operator<(FiltrationValue f1, FiltrationValue f2);
    /**
     * @brief Smaller or equal operator.
     */
    bool operator<=(FiltrationValue f1, FiltrationValue f2);
    /**
     * @brief Equality operator
     */
    bool operator==(FiltrationValue f1, FiltrationValue f2);
    /**
     * @brief Not equal operator
     */
    bool operator!=(FiltrationValue f1, FiltrationValue f2);

    //TODO: how expressing in the concept that a class should "overload" an external method?
    //Used the keyword `friend` to show that the class/method is not supposed to be inside FiltrationValue
    //but this can be misunderstood, as it does not have to be a friend.

    // overload of the standard numeric_limits class from #include <limits>
    friend class std::numeric_limits<FiltrationValue> {
     public:
      static constexpr bool has_infinity;

      //can simply throw when called if has_infinity == false
      static FiltrationValue infinity() throw();
      //can simply throw when called if has_infinity == true
      static FiltrationValue max() throw();

      //only used in decode_extended_filtration
      // static FiltrationValue quiet_NaN() throw();
    };

    /**
     * @brief Given two filtration values at which a simplex exists, returns the minimal union of births generating
     * a lifetime including those two values. The overload for native arithmetic types like `double` or `int` is
     * already implemented.
     * If the filtration is 1-critical and totally ordered (as in one parameter persistence), the union is simply
     * the minimum of the two values. If the filtration is 1-critical, but not totally ordered (possible for
     * multi-persistence), than the union is also the minium if the two given values are comparable and the 
     * method should throw an error if they are not, as a same simplex should not exist at those two values.
     * Finally, if the filtration is k-critical, FiltrationValue should be able to store an union of values and this
     * method adds the values of @p f2 in @p f1 and removes the values from @p f1 which are comparable and greater than
     * other values.
     */
    friend void unify_births(FiltrationValue& f1, const FiltrationValue& f2);

    /**
     * @brief Given two filtration values, stores in the first value the greatest common upper bound of the two values.
     * The overload for native arithmetic types like `double` or `int` is already implemented.
     */
    friend void push_to_greatest_common_upper_bound(FiltrationValue& f1, const FiltrationValue& f2);
  };
