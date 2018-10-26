package com.lucidlogic.horsetracker.adapter.viewholder;

import android.databinding.ViewDataBinding;
import android.support.v7.widget.RecyclerView;

public class GenericDataBindingViewHolder extends RecyclerView.ViewHolder {

    private ViewDataBinding mViewDataBinding;

    public GenericDataBindingViewHolder( ViewDataBinding viewDataBinding) {
        super(viewDataBinding.getRoot());

        mViewDataBinding = viewDataBinding;
        mViewDataBinding.executePendingBindings();
    }

    public ViewDataBinding getViewDataBinding() {
        return mViewDataBinding;
    }
}