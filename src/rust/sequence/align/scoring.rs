type Score = i32;

trait ScoringScheme<Symbol> {
    fn score(&self, symbol1: Symbol, symbol2: Symbol) -> Score;
}

struct SubstitutionMatrix<Symbol> {
    matrix: Vec<Vec<Score>>,
}

impl<Symbol> ScoringScheme<Symbol> for SubstitutionMatrix<Symbol> {
    fn score(&self, symbol1: Symbol, symbol2: Symbol) -> Score {
        self.matrix[symbol1][symbol2]
    }
}

struct Match<Symbol> {
    match_score: Score,
    mismatch_score: Score,
}

impl<Symbol> ScoringScheme<Symbol> for Match<Symbol> {
    fn score(&self, symbol1: Symbol, symbol2: Symbol) -> Score {
        if symbol1 == symbol2 {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}