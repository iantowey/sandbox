/** @file Doxy Test File */

struct mystruct1 {
  int a;
  int b;
};

struct mystruct2 {
  double x;
  double y;
};

/** @brief Do something to a mystruct1.
 *
 * What we do is a bit of a mystery, but trust
 * me, it works.
 *
 * @param s The struct mystruct1 to do it to.
 * @return an integer
 */
int do_something_to_struct1(struct mystruct1 *s);

/** @brief Do something to a mystruct2.
 *
 * Whatever it is, it sure is wonderful.
 *
 * @param s The struct mystruct2 to do it to.
 * @return a double
 */
double do_something_to_struct2(struct mystruct2 *s);
