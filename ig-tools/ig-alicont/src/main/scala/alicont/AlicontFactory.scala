package alicont

import alicont.algorithms.AlgorithmType
import alicont.algorithms.AlgorithmType.AlgorithmType
import alicont.conts.{AbstractAlicont, simple, affine}

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:05
 */
object AlicontFactory {

  def createSimpleAlicont(maxheight : Int, query : String, gap : Double, score_matrix : Array[Array[Double]],
                          algo_type : AlgorithmType) : AbstractAlicont = {
    algo_type match {
      case AlgorithmType.GLOBAL     => new simple.AlicontGlobal(maxheight, query, gap, score_matrix)
      case AlgorithmType.LOCAL      => new simple.AlicontLocal(maxheight, query, gap, score_matrix)
      case AlgorithmType.SEMIGLOBAL => new simple.AlicontSemiglobal(maxheight, query, gap, score_matrix)
      case _ => null
    }
  }

  def createAffineAlicont(maxheight : Int, query : String, gap_open : Double, gap_ext : Double,
                          score_matrix : Array[Array[Double]], algo_type : AlgorithmType) : AbstractAlicont = {
    algo_type match {
      case AlgorithmType.AFFINE_GLOBAL     => new affine.AlicontGlobal(maxheight, query, gap_open, gap_ext, score_matrix)
      case AlgorithmType.AFFINE_LOCAL      => new affine.AlicontLocal(maxheight, query, gap_open, gap_ext, score_matrix)
      case AlgorithmType.AFFINE_SEMIGLOBAL => new affine.AlicontSemiglobal(maxheight, query, gap_open, gap_ext, score_matrix)
      case _ => null
    }
  }
}
