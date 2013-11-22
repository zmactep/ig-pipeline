package common


/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 21.11.13
 * Time: 17:28
 */
object Algo {
  def longestSubstr(s : String, t : String) : (Int, Int, Int) = {
    if (s.isEmpty || t.isEmpty) {
      return (0, 0, 0)
    }

    val m = s.length
    val n = t.length
    var cost = 0
    var maxLen = 0
    var posi = 0
    var posj = 0
    var p = Array.fill[Int](n)(0)
    var d = Array.fill[Int](n)(0)

    for (i <- 0 until m) {
      for (j <- 0 until n) {
        // calculate cost/score
        if (s.charAt(i) != t.charAt(j)) {
          cost = 0
        } else {
          if ((i == 0) || (j == 0)) {
            cost = 1
          } else {
            cost = p(j - 1) + 1
          }
        }
        d(j) = cost

        if (cost > maxLen) {
          posi = i
          posj = j
          maxLen = cost
        }
      } // for {}

      val swap = p
      p = d
      d = swap
    }

    (posi - maxLen + 1, posj - maxLen + 1, maxLen)
  }

  def translateString(dna : String) : String = {
    val result = StringBuilder.newBuilder
    val codon_table = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
    val nucleo = "ACGT"
    dna.zipWithIndex.foreach(tpl => {
      val (_, i) = tpl
      if (i % 3 == 0 && dna.size - i > 2) {
        val pos = nucleo.indexOf(dna(i))*16 + nucleo.indexOf(dna(i+1))*4 + nucleo.indexOf(dna(i+2))
        result.append(codon_table(pos))
      }
    })
    result.toString()
  }

  def reverseComp(dna : String) : String = {
    val nuc = "ACGT"
    dna.reverseMap(ch => if (ch == 'N') 'N' else nuc(nuc.size - nuc.indexOf(ch) - 1))
  }
}
