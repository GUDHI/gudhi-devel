/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - YYYY/MM author: description of the modification.
 *         
 */

/** \brief The concept <CODE>FilteredMorseComplex</CODE> is a refinement of 
 * <CODE>FilteredComplex</CODE>, which describes the requirements for a type to 
 * implement a filtered cell complex with a Morse matching.
 */
struct FilteredMorseComplex : public FilteredComplex {
/** \brief Tells is a cell is critical.
 * 
 * \details Returns <CODE>true</CODE> iff the cell is critical. 
 */
  bool critical(Simplex_handle sh);
/** \brief Returns a handle to the cell an input cell <CODE>sh</CODE> is paired 
 * with in the Morse matching. Returns <CODE>sh</CODE> itself if the cell is 
 * critical.
 */
  Simplex_handle paired_with(Simplex_handle sh);
/** \brief Pair two cells together in a Morse matching. */
  void assign_pairing(Simplex_handle sh_t, Simplex_handle sh_s);
/** \brief Make a cell critical in a Morse matching. */
  void make_critical(Simplex_handle sh);
/** \brief Tells if two cells are paired together in a Morse matching.
 * 
 * \details <CODE>true</CODE> iff the two handles point to the two cells of a 
 * Morse pair, or if they point to a same cell that is critical. 
 */
  bool is_paired_with(Simplex_handle sh_t, Simplex_handle sh_s);
};