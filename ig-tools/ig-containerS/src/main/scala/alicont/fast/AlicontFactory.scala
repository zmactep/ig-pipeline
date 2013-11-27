package alicont.fast

import alicont.fast.algorithms.AlgorithmType
import alicont.fast.algorithms.AlgorithmType.AlgorithmType
import alicont.fast.conts.simple.{AlicontSemiglobal, AlicontLocal, AlicontGlobal}

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:05
 */
object AlicontFactory {

  def createSimpleAlicont(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]],
                          algo_type : AlgorithmType) : AbstractAlicont = {
    algo_type match {
      case AlgorithmType.GLOBAL     => new AlicontGlobal(maxheight, query, gap, score_matrix)
      case AlgorithmType.LOCAL      => new AlicontLocal(maxheight, query, gap, score_matrix)
      case AlgorithmType.SEMIGLOBAL => new AlicontSemiglobal(maxheight, query, gap, score_matrix)
      case _ => null
    }
  }

  def createAffineAlicont(maxheight : Int, query : String, gap_open : Int, gap_ext : Int,
                          score_matrix : Array[Array[Int]], algo_type : AlgorithmType) : AbstractAlicont = {
    algo_type match {
      case AlgorithmType.AFFINE_GLOBAL     => null
      case AlgorithmType.AFFINE_LOCAL      => null
      case AlgorithmType.AFFINE_SEMIGLOBAL => null
      case _ => null
    }
  }
}
