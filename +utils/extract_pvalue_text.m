function pvalue_text = extract_pvalue_text(pval, is_corrected)

if nargin<2
    is_corrected = 0;
end

if ~is_corrected
    if pval>=0.01
        pvalue_text = sprintf('p = %.2f', pval);
    elseif pval<0.01 && pval>=0.001
        pvalue_text = 'p < 0.01';
    elseif pval<0.001
        pvalue_text = 'p < 0.001';
    end
else
    if pval>=0.01
        pvalue_text = sprintf('p_{FWER} = %.2f', pval);
    elseif pval<0.01 && pval>=0.001
        pvalue_text = 'p_{FWER} < 0.01';
    elseif pval<0.001
        pvalue_text = 'p_{FWER} < 0.001';
    end
end