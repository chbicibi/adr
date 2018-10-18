module GA_Base
  use Base_Func

  implicit none

  public

  type, extends(TGroup) :: Group1
      
  end type Group1

  interface Crossover_Order
    procedure :: Crossover_Order_i
  end interface Crossover_Order

  interface Crossover_SinglePoint
    procedure :: Crossover_SinglePoint_i, Crossover_SinglePoint_r
  end interface Crossover_SinglePoint

  interface Crossover_Uniform
    procedure :: Crossover_Uniform_i, Crossover_Uniform_r
  end interface Crossover_Uniform

  interface Crossover_BLXalpha
    procedure :: Crossover_BLXalpha_r
  end interface Crossover_BLXalpha

  contains

  ! ----------------------------------------------------------------------------
  ! selection
  ! ----------------------------------------------------------------------------
  function SelectSingle_Roulette(score)
    integer                  :: SelectSingle_Roulette(2)
    type(THash), intent(in)  :: score(:)
    type(THash), allocatable :: score_roulette(:)
    integer                  :: size_score, n(2), i

    size_score = size(score)

    score_roulette        = score
    score_roulette%value  = [(score(i)%value * Random(), i = 1, size_score)]
    ! score_roulette%value  = [(-1d0 / (1d0 + EXP(5d0 * score(i)%value)) * Random(), i = 1, size_score)]
    score_roulette        = Sort(score_roulette)

    n                     = Shuffle(2)
    SelectSingle_Roulette = score_roulette(n)%key
  end function SelectSingle_Roulette

  function SelectSingle_Tournament(score)
    integer                 :: SelectSingle_Tournament(2)
    type(THash), intent(in) :: score(:)
    integer                 :: size_score, size_tournament, n(2)

    size_score      = size(score)
    size_tournament = size_score / 5

    n(1:1) = Sort(Shuffle(size_score, size_tournament))
    n(2:2) = Sort(Shuffle(size_score - 1, size_tournament))

    if (n(1) == n(2)) n(2) = n(2) + 1

    SelectSingle_Tournament = score(n)%key
  end function SelectSingle_Tournament

  function SelectPateto_t(num_indiv, rank_pareto)
    integer              :: SelectPateto_t(2)
    integer, intent(in)  :: num_indiv, rank_pareto(:)
    integer, allocatable :: list(:)
    integer              :: n(2), i

    if (count(rank_pareto == 2) > 2) then
      list = PACK([(i, i = 1, num_indiv)], rank_pareto == 2)
    else
      list = PACK([(i, i = 1, num_indiv)], rank_pareto == 2 .or. rank_pareto == 3)
    end if

    n              = Shuffle(size(list), 2)
    SelectPateto_t = list(n)
  end function SelectPateto_t

  function SelectPateto(num_indiv, rank_pareto)
    integer             :: SelectPateto(2)
    integer, intent(in) :: num_indiv, rank_pareto(:)
    integer             :: n(2), t, i

    i = 1
    do
      t = Random(num_indiv)
      if ((i == 1 .or. t /= n(1)) .and. rank_pareto(t) > 0 .and. Random() < 1d0 / (1d0 + EXP(rank_pareto(t) * 0.3d0))) then
        n(i) = t
        i = i + 1
        if (i > 2) exit
      end if
    end do

    SelectPateto = n
  end function SelectPateto

  ! ----------------------------------------------------------------------------
  ! crossover
  ! ----------------------------------------------------------------------------
  function Crossover_Order_i(chromos_in, num_cross) result (chromos_out)
    type(TIAry)              :: chromos_out(2)
    type(TIAry), intent(in)  :: chromos_in(:)
    integer,     intent(in)  :: num_cross
    integer,     allocatable :: dup(:), gene_cross(:), gene_cross2(:)
    integer                  :: num_chrom, num_dup, i

    num_chrom = size(chromos_in(1)%elm)
    dup       = PACK(chromos_in(1)%elm, [(any(chromos_in(1)%elm(i) == chromos_in(2)%elm), i = 1, num_chrom)])
    num_dup   = size(dup)

    chromos_out(1:2) = chromos_in

    if (num_dup <= 1 .or. num_cross <= 1) then
      return
    else if (num_dup <= num_cross) then
      gene_cross = dup
    else
      gene_cross = dup(Sort(Shuffle(num_dup, num_cross)))
      num_dup    = num_cross
    end if

    gene_cross2 = PACK(chromos_in(2)%elm, [(any(chromos_in(2)%elm(i) == gene_cross), i = 1, num_chrom)])

    do i = 1, num_dup
      where (chromos_in(1)%elm == gene_cross(i) ) chromos_out(1)%elm = gene_cross2(i)
      where (chromos_in(2)%elm == gene_cross2(i)) chromos_out(2)%elm = gene_cross(i)
    end do
  end function Crossover_Order_i

  function Crossover_SinglePoint_i(chromos_in) result (chromos_out)
    type(TIAry)              :: chromos_out(2)
    type(TIAry), intent(in)  :: chromos_in(:)
    integer                  :: num_chrom, p, i

    num_chrom   = size(chromos_in(1)%elm)
    p           = Random(num_chrom - 1)
    chromos_out = [(TIAry([chromos_in(i)%elm(1:p), chromos_in(3 - i)%elm(p + 1:num_chrom)]), i = 1, 2)]
  end function Crossover_SinglePoint_i

  function Crossover_SinglePoint_r(chromos_in) result (chromos_out)
    type(TRAry)             :: chromos_out(2)
    type(TRAry), intent(in) :: chromos_in(:)
    integer                 :: num_chrom, p, i

    num_chrom   = size(chromos_in(1)%elm)
    p           = Random(num_chrom - 1)
    chromos_out = [(TRAry([chromos_in(i)%elm(1:p), chromos_in(3 - i)%elm(p + 1:num_chrom)]), i = 1, 2)]
  end function Crossover_SinglePoint_r

  function Crossover_Uniform_i(chromos_in) result (chromos_out)
    type(TIAry)              :: chromos_out(2)
    type(TIAry), intent(in)  :: chromos_in(:)
    integer,     allocatable :: locus(:)
    integer                  :: num_parent, num_chrom, i, j

    num_parent  = size(chromos_in)
    num_chrom   = size(chromos_in(1)%elm)
    locus       = [(Random(num_parent), i = 1, num_chrom)] - 1
    chromos_out = [(TIAry([(chromos_in(mod(locus(i) + j, num_parent) + 1)%elm(i), i = 1, num_chrom)]), j = 0, 1)]
  end function Crossover_Uniform_i

  function Crossover_Uniform_r(chromos_in) result (chromos_out)
    type(TRAry)              :: chromos_out(2)
    type(TRAry), intent(in)  :: chromos_in(:)
    integer,     allocatable :: locus(:)
    integer                  :: num_parent, num_chrom, i, j

    num_parent  = size(chromos_in)
    num_chrom   = size(chromos_in(1)%elm)
    locus       = [(Random(num_parent), i = 1, num_chrom)] - 1
    chromos_out = [(TRAry([(chromos_in(mod(locus(i) + j, num_parent) + 1)%elm(i), i = 1, num_chrom)]), j = 0, 1)]
  end function Crossover_Uniform_r

  function Crossover_BLXalpha_r(chromos_in, alpha) result (chromos_out)
    type(TRAry2)              :: chromos_out(1)
    type(TRAry2), intent(in)  :: chromos_in(:)
    real(8),      intent(in)  :: alpha
    real(8),      allocatable :: d(:, :), ave(:, :)

    d                          = abs(chromos_in(1)%elm - chromos_in(2)%elm) * (1d0 + 2d0 * alpha)
    ave                        = (chromos_in(1)%elm + chromos_in(2)%elm) / 2d0
    chromos_out(1)%elm         =     max(ave + Random(d, 0d0), 0d0)
    chromos_out(1)%elm(1, :)   =     max(chromos_out(1)%elm(1, :), 0.0d0)
    chromos_out(1)%elm(2, :)   =     max(chromos_out(1)%elm(2, :), 0.8d0)
    chromos_out(1)%elm(3, :)   = min(max(chromos_out(1)%elm(3, :), 0.1d0), 0.9d0)
    chromos_out(1)%elm(4:5, :) = min(max(chromos_out(1)%elm(4:5, :), 0d0), 1d0)
    chromos_out(1)%elm(6:7, :) =     max(chromos_out(1)%elm(6:7, :), 0d0)
  end function Crossover_BLXalpha_r

  ! ----------------------------------------------------------------------------
  ! mutation
  ! ----------------------------------------------------------------------------
  function Mutate_Swap(chrom) result (chrom_out)
    integer, allocatable :: chrom_out(:)
    integer, intent(in)  :: chrom(:)
    integer              :: num_chrom, locus(2)

    num_chrom = size(chrom)
    chrom_out = chrom
    locus     = Shuffle(num_chrom, 2)

    chrom_out(locus(1)) = chrom(locus(2))
    chrom_out(locus(2)) = chrom(locus(1))
  end function Mutate_Swap

  function Mutate_Change(chrom, num_target_all) result (chrom_out)
    integer, allocatable :: chrom_out(:)
    integer, intent(in)  :: chrom(:), num_target_all
    integer              :: num_chrom, locus_from, locus_to, p

    num_chrom  = size(chrom)
    chrom_out  = chrom

    p          = Random(num_chrom)
    locus_from = chrom(p)
    locus_to   = Random(num_target_all - 1)

    if (locus_to == locus_from) locus_to = locus_to + 1

    chrom_out(p) = 0
    where (chrom_out == locus_to) chrom_out = locus_from
    where (chrom_out == 0)        chrom_out = locus_to
  end function Mutate_Change

  function Mutate_Neighbor(chrom) result (chrom_out)
    real(8), allocatable :: chrom_out(:, :)
    real(8), intent(in)  :: chrom(:, :)
    integer              :: p, q

    chrom_out  = chrom
    p          = Random(ubound(chrom, 1))
    q          = Random(ubound(chrom, 2))

    chrom_out(p, q) = chrom_out(p, q) * Random(0.1d0, 1d0)

    select case (p)
    case (1, 6:7)
      chrom_out(p, q) = max(chrom_out(p, q), 0d0)
    case (2)
      chrom_out(p, q) = max(chrom_out(p, q), 0.8d0)
    case (3)
      chrom_out(p, q) = min(max(chrom_out(p, q), 0.1d0), 0.9d0)
    case (4:5)
      chrom_out(p, q) = min(max(chrom_out(p, q), 0d0), 1d0)
      chrom_out(p, q) = max(chrom_out(p, q), 0d0)
      chrom_out(p, q) = max(chrom_out(p, q), 0d0)
    end select
  end function Mutate_Neighbor
end module GA_Base
