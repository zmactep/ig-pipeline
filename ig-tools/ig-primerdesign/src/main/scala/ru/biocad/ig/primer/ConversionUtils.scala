package ru.biocad.ig.primer

import DnaUtils.Sequence
import Predef.{augmentString => _, _}

/**
 * Created by Kos on 03.03.14.
 */
/**
 * See http://stackoverflow.com/questions/1886953/is-there-a-way-to-control-which-implicit-conversion-will-be-the-default-used/1887678#1887678
 * to understand "extends LowPriorityImplicits": "if alternatives have identical argument types, the one which is defined in a subclass wins".
 * Otherwise LowPriorityImplicits has wrapString() implicit conversion, which is ambiguous with String2Sequence()
 */
object ConversionUtils extends LowPriorityImplicits {
  /**
   * implicit conversion from sequence string to Sequence (Used in DnaUtils).
   * @param sequence input sequence like "ACT"
   * @return Sequence like List(Set("A"), Set("C"), Set("T"))
   */
  implicit def String2Sequence(sequence: String): Sequence = sequence.toCharArray.map(c => Set(c.toString)).toList

  /**
   * Adds implicit possibility to modify i-th Set[String] in List
   * @param s Sequence
   */
  implicit class SequenceImprovements(s: Sequence) {
    def replaceAt(i: Int, replacement: Set[String]): Sequence = s.updated(i, replacement)
    def addAlternativeAt(i: Int, alternative: Set[String]): Sequence = s.replaceAt(i, alternative ++ s(i))
  }
}

