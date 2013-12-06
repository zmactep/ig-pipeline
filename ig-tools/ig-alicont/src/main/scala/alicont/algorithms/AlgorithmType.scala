package alicont.algorithms

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:16
 */
object AlgorithmType extends Enumeration {
  type AlgorithmType = Value
  val GLOBAL, LOCAL, SEMIGLOBAL = Value
  val AFFINE_GLOBAL, AFFINE_LOCAL, AFFINE_SEMIGLOBAL = Value

  def affine : List[AlgorithmType] = AFFINE_GLOBAL :: AFFINE_LOCAL :: AFFINE_SEMIGLOBAL :: Nil

  def simple : List[AlgorithmType] = GLOBAL :: LOCAL :: SEMIGLOBAL :: Nil
}
